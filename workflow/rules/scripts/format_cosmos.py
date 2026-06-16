"""Format per-condition variant scores into a single cosmos input CSV.

cosmos (pimentellab/cosmos, `pip install cosmos-dms`) models direct vs.
indirect effects across two (or more) sequential phenotypes. Its `DMSData`
loader requires a single wide table with, per `DMSData._check_cols`
(src/cosmos/dms_data/dms_data.py):

    {"variants", "group", "type"}
    | {f"beta_hat_{i+1}" for i in range(len(phenotypes))}
    | {f"se_hat_{i+1}"   for i in range(len(phenotypes))}

The numeric `_N` suffix is the ENTIRE phenotype interface — the phenotype
*names* the analyst later passes to `DMSData([pheno1, pheno2])` are used only
for their count and as plot labels, never matched against column names. So
dumpling's job is purely: assign each experimental condition a stable slot
number N (declared via the optional `phenotype` column in the experiment CSV)
and emit its score -> `beta_hat_N`, sd -> `se_hat_N`.

Column mapping from a dumpling score CSV (rosace example: variants, position,
wildtype, mutation, type, mean, sd, ...):

    variants  -> variants   (join key across conditions; emitted as-is)
    position  -> group      (cosmos sorts + bins this; must be sortable)
    type      -> type       (combined with mutation length, see derive_cosmos_type)
    mean      -> beta_hat_N
    sd        -> se_hat_N

Join policy: OUTER join on `variants`. A variant scored in one phenotype but
not another keeps its row with blank beta_hat/se_hat for the missing
phenotype(s). This mirrors format_mavedb's "surface filter outcomes explicitly,
don't silently lose signal" convention. `DMSData._check_cols` only checks
column presence, so the NA-padded file loads fine; whether downstream prior
fitting tolerates NA rows is a separate, open question (see
docs/cosmos_export_design.md).

Run via Snakemake or standalone:
  python format_cosmos.py --scores cond_A.csv cond_B.csv \\
      --output exp_cosmos.csv [--backend rosace]
  # --scores are given in slot order: first -> beta_hat_1, second -> beta_hat_2.
"""

import argparse

import pandas as pd

# Default column names per scoring backend. Mirrors format_mavedb._BACKEND_DEFAULTS
# but adds the columns cosmos needs that MaveDB didn't (position->group, and the
# type/mutation pair that drives derive_cosmos_type).
#
# Rosace verified against results/example_experiment/rosace/cond_A_scores.csv:
#   variants, position, wildtype, mutation, type, mean, sd, lfsr, ...
# Lilace: variant, type, <metadata>, position, effect, effect_se, lfsr, ...
#   NOTE: lilace's per-variant indel length column is unverified — confirm it
#   emits an `I_N`/`D_N`-style `mutation` column before trusting derive_cosmos_type
#   on lilace indels. (See docs/cosmos_export_design.md open questions.)
_BACKEND_DEFAULTS = {
    "rosace": {
        "variants": "variants",
        "score": "mean",
        "sd": "sd",
        "position": "position",
        "type": "type",
        "mutation": "mutation",
    },
    "lilace": {
        "variants": "variant",
        "score": "effect",
        "sd": "effect_se",
        "position": "position",
        "type": "type",
        "mutation": "mutation",
    },
}


def derive_cosmos_type(type_label: str, mutation: str) -> str:
    """Map a dumpling (type, mutation) pair to a single cosmos `type` label.

    cosmos matches `type` against the analyst's `include_type` / `exclude_type`
    lists as plain strings — there is NO fixed vocabulary and NO length cap, so
    you are free to emit whatever granularity your library actually contains.

    Inputs (from the dumpling score CSV):
      type_label : one of "missense", "synonymous", "nonsense", "insertion",
                   "deletion", "wt"  (the `type` column)
      mutation   : for substitutions, the substituted AA; for indels, a
                   length-encoded code from process_variants.py —
                     "I_<len>"  insertion of <len> codons   (process_variants.py:105)
                     "D_<len>"  deletion of <len> codons     (:223)
                     "ID_<len>" in-frame insertion-deletion  (:248)

    This is the one transform in the whole exporter that encodes a real
    decision rather than a mechanical rename, which is why it's left to you:

      1. Spelling of length-resolved indels. cosmos's vignette uses
         "insertion3" / "deletion2" (no separator). Does that match the
         `include_type`/`exclude_type` lists you'll pass to DMSData, or do you
         prefer "insertion_3"? It only has to agree with your own analysis code.
      2. In-frame insdels ("ID_<len>"): a distinct cosmos type, fold into
         insertion/deletion, or exclude entirely?
      3. Pass-throughs: presumably missense/synonymous/nonsense/wt map to
         themselves unchanged — confirm that's what your include/exclude lists
         expect.

    Return the cosmos `type` string for this variant.

    Decisions made here (the open choices from the docstring above):
      1. Spelling: no separator ("insertion3"), matching cosmos's own vignette.
      2. In-frame insdels ("ID_<len>") keep a distinct, length-resolved label
         (the variant's `type` string + length) rather than being folded into
         insertion/deletion — most informative, and the analyst can include or
         exclude it explicitly.
      3. Indels are detected by the authoritative mutation code prefix
         (I_/D_/ID_ from process_variants.py), not the `type` string, so a
         backend that labels its `type` column differently doesn't break this.
      4. Substitutions (missense/synonymous/nonsense) and wt pass through.
    """
    mutation = "" if mutation is None or mutation != mutation else str(mutation)
    head = mutation.split("_")[0] if "_" in mutation else ""
    if head in ("I", "D", "ID"):
        # Fuse the variant's type label with the codon length, no separator.
        return f"{type_label}{mutation.split('_')[-1]}"
    # Substitutions and wt: the type label is already the cosmos label.
    return type_label


def format_cosmos(score_paths, backend, cols, output_path):
    """Build the wide cosmos CSV from per-condition score CSVs.

    score_paths : list of paths in SLOT ORDER. score_paths[i] becomes
                  beta_hat_{i+1} / se_hat_{i+1}.
    cols        : dict of source column names (see _BACKEND_DEFAULTS).
    """
    beta_frames = []
    # group/type are properties of the variant, not the condition, so we
    # coalesce them across conditions: take the first non-null we see. On an
    # outer join a variant missing from condition 1 still gets its group/type
    # from whichever condition did score it.
    meta = None

    for i, path in enumerate(score_paths):
        n = i + 1
        df = pd.read_csv(path)

        required = [
            cols["variants"],
            cols["score"],
            cols["sd"],
            cols["position"],
            cols["type"],
            cols["mutation"],
        ]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                f"Expected columns not found in {path}: {missing}\n"
                f"Available columns: {list(df.columns)}\n"
                f"Override defaults via the backend column config."
            )

        key = df[cols["variants"]]
        if key.duplicated().any():
            dupes = key[key.duplicated()].unique().tolist()
            raise ValueError(
                f"{path}: duplicate variants in a single condition: {dupes[:5]}"
                f"{' ...' if len(dupes) > 5 else ''}. cosmos needs one row per "
                "variant per phenotype."
            )

        beta_frames.append(
            pd.DataFrame(
                {
                    "variants": key,
                    f"beta_hat_{n}": df[cols["score"]],
                    f"se_hat_{n}": df[cols["sd"]],
                }
            ).set_index("variants")
        )

        this_meta = pd.DataFrame(
            {
                "variants": key,
                "group": df[cols["position"]],
                "type": [
                    derive_cosmos_type(t, m)
                    for t, m in zip(df[cols["type"]], df[cols["mutation"]])
                ],
            }
        ).set_index("variants")
        meta = this_meta if meta is None else meta.combine_first(this_meta)

    # axis=1 concat aligns on the index and defaults to an OUTER join -> the
    # union of variants across all conditions, NA-padded where a variant is
    # missing from a phenotype. This IS the decided join policy.
    wide = pd.concat(beta_frames, axis=1)
    wide = meta.join(wide, how="right")  # attach group/type to every scored row

    # Column order: variants, group, type, then beta_hat_1/se_hat_1, ...
    ordered = ["group", "type"]
    for i in range(len(score_paths)):
        ordered += [f"beta_hat_{i+1}", f"se_hat_{i+1}"]
    wide = wide[ordered]

    wide.reset_index().rename(columns={"index": "variants"}).to_csv(
        output_path, index=False
    )


def main():
    # Snakemake injects `snakemake` as a module global (not visible to dir());
    # check globals() — same gotcha documented in format_mavedb.main().
    if "snakemake" in globals():
        backend = snakemake.params.backend
        defaults = _BACKEND_DEFAULTS.get(backend, _BACKEND_DEFAULTS["rosace"])
        cfg = snakemake.config.get("cosmos", {})
        cols = {k: cfg.get(f"{k}_column", v) for k, v in defaults.items()}
        format_cosmos(
            score_paths=list(snakemake.input.scores),
            backend=backend,
            cols=cols,
            output_path=snakemake.output.cosmos,
        )
        return

    parser = argparse.ArgumentParser(description="Format scores for cosmos.")
    parser.add_argument(
        "--scores",
        nargs="+",
        required=True,
        help="Per-condition score CSVs in slot order (1st -> beta_hat_1, ...).",
    )
    parser.add_argument("--output", required=True, help="Output cosmos CSV path")
    parser.add_argument("--backend", default="rosace", choices=list(_BACKEND_DEFAULTS))
    args = parser.parse_args()

    format_cosmos(
        score_paths=args.scores,
        backend=args.backend,
        cols=_BACKEND_DEFAULTS[args.backend],
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
