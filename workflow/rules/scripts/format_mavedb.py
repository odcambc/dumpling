"""Format variant scores and raw counts into MaveDB score/count CSV format.

Score output (`{cond}_mavedb.csv`): columns hgvs_pro, score, sd (if available).
Count output (`{cond}_mavedb_counts.csv`): columns hgvs_pro plus one column per
input sample (sample name as-is from the dumpling experiment CSV).

The hgvs_pro column is emitted in MAVE-HGVS form — three-letter amino acid
codes, no parenthesized "predicted" wrapper, `Ter` for stop. Rosace and
lilace emit one-letter parenthesized form (`p.(M1L)`, `p.(A1*)`);
to_mave_hgvs handles the structural conversion. Verified against the
MAVE-HGVS spec at github.com/VariantEffect/mavehgvs/docs/spec.rst
(audit 2026-05-26):

  rosace/lilace          MAVE-HGVS
  -----------------      ----------------
  p.(M1L)                p.Met1Leu
  p.(M1=)                p.Met1=
  p.(A1del)              p.Ala1del
  p.(A1_C3del)           p.Ala1_Cys3del
  p.(A1_A2insGGG)        p.Ala1_Ala2insGlyGlyGly
  p.(A1*)                p.Ala1Ter

Score and count tables MUST share the same variant set per the MaveDB spec
("Score and count tables must share the same set of variants and the same
index column"). We build the union of variants across the score CSV and all
per-sample count CSVs for the condition, then emit both files indexed on
that union. Variants observed in counts but filtered out by rosace/lilace
appear in the score CSV with NaN score/sd — depositors see the filter
outcomes explicitly rather than silently losing the raw signal.

Run via Snakemake or standalone:
  python format_mavedb.py --scores <scores.csv> --score-output <out.csv> \\
      [--counts <a.csv> <b.csv> ...] [--sample-names A B ...] \\
      [--count-output <counts.csv>] [--backend rosace]
"""

import argparse
import re

import pandas as pd

# Optional: mavehgvs validates the converted string is in MAVE-HGVS form.
# Available in the dumpling_env conda env (added via pip). If not importable
# (e.g. user ran without --use-conda), we skip validation rather than break
# the deposit path — the structural conversion is still applied.
try:
    from mavehgvs import Variant as _MaveVariant  # type: ignore
except ImportError:
    _MaveVariant = None

# Default column names per scoring backend. Override via config:
# mavedb.hgvs_column, mavedb.score_column, mavedb.sd_column.
#
# Rosace verified against results/example_experiment/rosace/cond_A_scores.csv:
#   variants, position, wildtype, mutation, type, mean, sd, lfsr, ...
# Lilace verified against pimentellab/lilace R/lilace.R `.get_score_df()`:
#   variant, type, <metadata>, position, effect, effect_se, lfsr, ...
# Note the schema drift: rosace uses plural `variants` + `mean`/`sd`; lilace
# uses singular `variant` + `effect`/`effect_se`. Don't merge these into a
# shared default — the column names ARE different.
_BACKEND_DEFAULTS = {
    "rosace": {"hgvs_col": "variants", "score_col": "mean", "sd_col": "sd"},
    "lilace": {"hgvs_col": "variant", "score_col": "effect", "sd_col": "effect_se"},
}

_SYNONYMOUS_RE = re.compile(r"^p\.\(([A-Z])(\d+)([A-Z])\)$")

# One-letter to three-letter amino acid mapping for MAVE-HGVS conversion.
# Inverted from process_variants.aa_3to1_dict. Both `X` (dumpling's
# internal stop code from aa_3to1_dict["Stp"] -> "X", which is what shows
# up in rosace/lilace score CSVs — verified empirically against
# results/example_experiment/rosace/cond_A_scores.csv 2026-05-27) and `*`
# (standard one-letter HGVS stop notation) map to "Ter", the only stop
# code MAVE-HGVS accepts. The spec explicitly rejects `*` and `X` literal
# stops in protein variants; the canonical form is always `Ter`.
_AA_1to3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "X": "Ter",  # dumpling's internal stop code (process_variants.aa_3to1_dict)
    "*": "Ter",  # alternate stop notation, kept for robustness
}

# Inner-content patterns (after stripping the outer `p.( ... )`).
_SUBSTITUTION = re.compile(r"^([A-Z*])(\d+)([A-Z*=])$")
_SINGLE_DELETION = re.compile(r"^([A-Z])(\d+)del$")
_MULTI_DELETION = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)del$")
_INSERTION = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")
_STRIP_OUTER_PARENS = re.compile(r"^p\.\((.+)\)$")


def normalize_hgvs(hgvs: str) -> str:
    """Convert same-AA synonymous notation p.(X1X) → p.(X1=) per HGVS standard."""
    m = _SYNONYMOUS_RE.match(hgvs)
    if m and m.group(1) == m.group(3):
        return f"p.({m.group(1)}{m.group(2)}=)"
    return hgvs


def _expand_aa_seq(seq: str) -> str:
    """Convert a one-letter AA sequence (`GGG`) to MAVE-HGVS three-letter
    concatenation (`GlyGlyGly`). Raises on unknown codes."""
    out = []
    for ch in seq:
        if ch not in _AA_1to3:
            raise ValueError(
                f"Unknown amino acid one-letter code {ch!r} in {seq!r}"
            )
        out.append(_AA_1to3[ch])
    return "".join(out)


def to_mave_hgvs(hgvs: str) -> str:
    """Convert rosace/lilace HGVS (`p.(M1L)`, `p.(A1_C3del)`, etc.) to the
    MAVE-HGVS canonical form (three-letter AA codes, no parenthesized
    "predicted" wrapper, `Ter` for stop).

    MAVE-HGVS rejects `p.(...)` parenthesization for everything except the
    population-level synonymous form `p.(=)`. Rosace produces one-letter
    parenthesized output for ALL variant classes; this function does the
    structural conversion. Returns the input unchanged when it's already
    `p.(=)`.

    Raises ValueError for shapes we don't recognize — better to fail-loud at
    format time than ship a MaveDB CSV with rows that'll be rejected at
    deposit. Catches the upstream-issue #31 mode where users discovered
    the original deposit produced unusable files.
    """
    if hgvs == "p.(=)":
        return hgvs  # Population-level synonymous; the one valid p.(...) form.

    m = _STRIP_OUTER_PARENS.match(hgvs)
    inner = m.group(1) if m else (hgvs[2:] if hgvs.startswith("p.") else hgvs)

    # Substitution / synonymous-at-position / nonsense.
    m = _SUBSTITUTION.match(inner)
    if m:
        wt, pos, alt = m.group(1), m.group(2), m.group(3)
        wt3 = _AA_1to3[wt]
        if alt == "=":
            return f"p.{wt3}{pos}="
        return f"p.{wt3}{pos}{_AA_1to3[alt]}"

    # Single-residue deletion.
    m = _SINGLE_DELETION.match(inner)
    if m:
        return f"p.{_AA_1to3[m.group(1)]}{m.group(2)}del"

    # Multi-residue deletion.
    m = _MULTI_DELETION.match(inner)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), m.group(2), m.group(3), m.group(4)
        return f"p.{_AA_1to3[wt1]}{pos1}_{_AA_1to3[wt2]}{pos2}del"

    # Insertion (the inserted residues are themselves one-letter coded).
    m = _INSERTION.match(inner)
    if m:
        wt1, pos1, wt2, pos2, inserted = m.groups()
        return (
            f"p.{_AA_1to3[wt1]}{pos1}_{_AA_1to3[wt2]}{pos2}"
            f"ins{_expand_aa_seq(inserted)}"
        )

    raise ValueError(
        f"Unrecognized HGVS shape for MAVE-HGVS conversion: {hgvs!r}. "
        "Expected one of: missense (p.(A1V)), synonymous (p.(A1=)), "
        "single deletion (p.(A1del)), multi deletion (p.(A1_C3del)), "
        "insertion (p.(A1_A2insGGG)), or nonsense (p.(A1*))."
    )


def _validate_mave_hgvs(hgvs: str, row_idx: int) -> str:
    """If mavehgvs is importable, validate the string. Re-raise with row
    context on failure so users can find the offending row in the source CSV.
    No-op when mavehgvs isn't installed (returns the string unchanged)."""
    if _MaveVariant is None:
        return hgvs
    try:
        _MaveVariant(hgvs)
    except Exception as exc:  # noqa: BLE001 - mavehgvs raises a custom class
        raise ValueError(
            f"Row {row_idx}: MAVE-HGVS validation failed for {hgvs!r}: {exc}. "
            "If mavehgvs rejected output our to_mave_hgvs converter produced, "
            "the converter's coverage of HGVS shapes needs updating."
        ) from exc
    return hgvs


def _to_mave_key_series(hgvs_series):
    """Apply normalize_hgvs + to_mave_hgvs to a series of dumpling-style HGVS
    strings, returning a series of MAVE-HGVS keys. Shared by scores and counts
    paths so the same input string always maps to the same index value."""
    return hgvs_series.apply(normalize_hgvs).apply(to_mave_hgvs)


def _aggregate_sample_counts(counts_path, count_col="count", hgvs_col="hgvs"):
    """Read one processed_counts/{sample}.csv and return a Series mapping
    MAVE-HGVS variant -> summed count.

    dumpling emits one row per (variant, codon) in processed_counts. Multiple
    codons can encode the same protein variant (e.g. M1F from both TTC and
    TTT), so we sum across the protein-level key. Synonymous codons all
    collapse to the same p.{wt}{pos}= key after normalize_hgvs."""
    df = pd.read_csv(counts_path)
    if count_col not in df.columns or hgvs_col not in df.columns:
        raise ValueError(
            f"{counts_path}: expected columns {count_col!r} and {hgvs_col!r}; "
            f"got {list(df.columns)}"
        )
    keys = _to_mave_key_series(df[hgvs_col])
    return df[count_col].groupby(keys).sum()


def format_scores(
    scores_path,
    backend,
    hgvs_col,
    score_col,
    sd_col,
    output_path,
    counts_paths=None,
    sample_names=None,
    counts_output_path=None,
):
    """Emit MaveDB score CSV, and optionally a matching count CSV.

    When counts_paths is provided, builds the union of variants across the
    score CSV and all per-sample count CSVs, then emits both files indexed on
    that union (NaN score/sd for variants only present in counts). When
    counts_paths is None, behaves as the original scores-only formatter.
    """
    scores = pd.read_csv(scores_path)

    missing = [c for c in [hgvs_col, score_col] if c not in scores.columns]
    if missing:
        raise ValueError(
            f"Expected columns not found in {scores_path}: {missing}\n"
            f"Available columns: {list(scores.columns)}\n"
            f"Set mavedb.hgvs_column / mavedb.score_column in config to override defaults."
        )

    # Two-step conversion: first canonicalize rosace's same-AA missense to
    # synonymous notation (p.(A1A) -> p.(A1=)), then convert the whole row
    # to MAVE-HGVS form (three-letter codes, no outer parens, Ter for stop).
    score_keys = _to_mave_key_series(scores[hgvs_col])
    scored_df = pd.DataFrame({"hgvs_pro": score_keys, "score": scores[score_col]})
    if sd_col and sd_col in scores.columns:
        scored_df["sd"] = scores[sd_col]
    scored_df = scored_df.set_index("hgvs_pro")

    emit_counts = counts_paths and counts_output_path
    if emit_counts:
        if not sample_names or len(sample_names) != len(counts_paths):
            raise ValueError(
                "sample_names must align 1:1 with counts_paths "
                f"(got {len(sample_names) if sample_names else 0} names "
                f"for {len(counts_paths)} files)"
            )

        # Per-sample MAVE-HGVS -> count series, joined column-wise.
        per_sample = {
            name: _aggregate_sample_counts(path)
            for name, path in zip(sample_names, counts_paths)
        }
        counts_df = pd.DataFrame(per_sample).fillna(0).astype(int)
        # Preserve the sample order Snakemake passed in (DataFrame() may
        # reorder by insertion-iteration order from the dict above, which is
        # insertion order in Py3.7+, but being explicit is safer).
        counts_df = counts_df[list(sample_names)]

        union_index = scored_df.index.union(counts_df.index)
        scored_df = scored_df.reindex(union_index)  # NaN-pads dropped rows
        counts_df = counts_df.reindex(union_index, fill_value=0)
    else:
        union_index = scored_df.index

    # Validate every emitted key (cheap; runs once per row).
    for i, key in enumerate(union_index):
        _validate_mave_hgvs(key, i)

    # NaN-padded rows (variants present only in counts, dropped by rosace's
    # filter) serialize as empty cells under pandas' default to_csv behavior:
    #     p.Ala2Val,,
    # The MaveDB spec doesn't document whether empty score/sd cells are
    # accepted at deposit validation. If upload rejects them, switch to
    # `to_csv(output_path, index=False, na_rep="NA")` so the empties become
    # explicit "NA" strings — or drop the NaN rows entirely from the score
    # CSV, at the cost of variant-set misalignment with the counts CSV
    # (which would re-trigger the spec violation this code exists to avoid).
    scored_df.reset_index().rename(columns={"index": "hgvs_pro"}).to_csv(
        output_path, index=False
    )
    if emit_counts:
        counts_df.reset_index().rename(columns={"index": "hgvs_pro"}).to_csv(
            counts_output_path, index=False
        )


def main():
    # Snakemake execution path. Snakemake's `script:` directive injects
    # `snakemake` as a module-level global. `dir()` inside a function
    # returns LOCAL names, not module globals — so `"snakemake" in dir()`
    # always returns False here and the original guard silently fell
    # through to the argparse path. Check globals() instead. (Found
    # 2026-05-27 when wiring deposit_to_mavedb default-on into get_input
    # and running the rule via snakemake.)
    if "snakemake" in globals():
        cfg = snakemake.config.get("mavedb", {})
        backend = snakemake.params.backend
        defaults = _BACKEND_DEFAULTS.get(backend, _BACKEND_DEFAULTS["rosace"])

        # Counts inputs are optional: when the rule declares them they appear
        # as snakemake.input.counts (list) and snakemake.params.sample_names
        # (list, aligned by index). Older rule wiring without counts still
        # works — format_scores falls back to scores-only.
        counts_paths = list(getattr(snakemake.input, "counts", []) or [])
        sample_names = list(getattr(snakemake.params, "sample_names", []) or [])
        counts_output_path = (
            snakemake.output.counts if counts_paths else None
        ) if hasattr(snakemake.output, "counts") else None

        format_scores(
            scores_path=snakemake.input.scores,
            backend=backend,
            hgvs_col=cfg.get("hgvs_column", defaults["hgvs_col"]),
            score_col=cfg.get("score_column", defaults["score_col"]),
            sd_col=cfg.get("sd_column", defaults["sd_col"]),
            output_path=snakemake.output.scores,
            counts_paths=counts_paths,
            sample_names=sample_names,
            counts_output_path=counts_output_path,
        )
        return

    # Standalone CLI path
    parser = argparse.ArgumentParser(description="Format scores/counts for MaveDB.")
    parser.add_argument("--scores", required=True, help="Score CSV from rosace or lilace")
    parser.add_argument("--score-output", required=True, help="Output MaveDB score CSV path")
    parser.add_argument("--backend", default="rosace", choices=list(_BACKEND_DEFAULTS))
    parser.add_argument("--hgvs-column", default=None)
    parser.add_argument("--score-column", default=None)
    parser.add_argument("--sd-column", default=None)
    parser.add_argument(
        "--counts", nargs="*", default=[],
        help="Per-sample processed_counts CSVs for the condition (optional).",
    )
    parser.add_argument(
        "--sample-names", nargs="*", default=[],
        help="Sample names aligned 1:1 with --counts. Used as count column headers.",
    )
    parser.add_argument(
        "--count-output", default=None,
        help="Output MaveDB count CSV path. Required when --counts is given.",
    )
    args = parser.parse_args()

    if args.counts and not args.count_output:
        parser.error("--count-output is required when --counts is given.")

    defaults = _BACKEND_DEFAULTS[args.backend]
    format_scores(
        scores_path=args.scores,
        backend=args.backend,
        hgvs_col=args.hgvs_column or defaults["hgvs_col"],
        score_col=args.score_column or defaults["score_col"],
        sd_col=args.sd_column or defaults["sd_col"],
        output_path=args.score_output,
        counts_paths=args.counts or None,
        sample_names=args.sample_names or None,
        counts_output_path=args.count_output,
    )


if __name__ == "__main__":
    main()
