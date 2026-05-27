"""Format variant scores into MaveDB score set CSV format.

Output columns: hgvs_pro, score, sd (if available). The hgvs_pro column is
emitted in MAVE-HGVS form — three-letter amino acid codes, no parenthesized
"predicted" wrapper, `Ter` for stop. Rosace and lilace emit one-letter
parenthesized form (`p.(M1L)`, `p.(A1*)`); to_mave_hgvs handles the
structural conversion. Verified against the MAVE-HGVS spec at
github.com/VariantEffect/mavehgvs/docs/spec.rst (audit 2026-05-26):

  rosace/lilace          MAVE-HGVS
  -----------------      ----------------
  p.(Met1Leu)            (n/a — already three-letter? no — rosace is 1-letter)
  p.(M1L)                p.Met1Leu
  p.(M1=)                p.Met1=
  p.(A1del)              p.Ala1del
  p.(A1_C3del)           p.Ala1_Cys3del
  p.(A1_A2insGGG)        p.Ala1_Ala2insGlyGlyGly
  p.(A1*)                p.Ala1Ter

Run via Snakemake or standalone:
  python format_mavedb.py --scores <scores.csv> --output <out.csv> [--backend rosace]
"""

import argparse
import re
import sys

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


def format_scores(scores_path, backend, hgvs_col, score_col, sd_col, output_path):
    scores = pd.read_csv(scores_path)

    missing = [c for c in [hgvs_col, score_col] if c not in scores.columns]
    if missing:
        raise ValueError(
            f"Expected columns not found in {scores_path}: {missing}\n"
            f"Available columns: {list(scores.columns)}\n"
            f"Set mavedb.hgvs_column / mavedb.score_column in config to override defaults."
        )

    out = pd.DataFrame()
    # Two-step conversion: first canonicalize rosace's same-AA missense to
    # synonymous notation (p.(A1A) -> p.(A1=)), then convert the whole row
    # to MAVE-HGVS form (three-letter codes, no outer parens, Ter for stop).
    # Validation runs per-row so error messages can name the offending index.
    normalized = scores[hgvs_col].apply(normalize_hgvs).apply(to_mave_hgvs)
    out["hgvs_pro"] = [
        _validate_mave_hgvs(s, i) for i, s in enumerate(normalized)
    ]
    out["score"] = scores[score_col]

    if sd_col and sd_col in scores.columns:
        out["sd"] = scores[sd_col]

    out.to_csv(output_path, index=False)


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

        format_scores(
            scores_path=snakemake.input.scores,
            backend=backend,
            hgvs_col=cfg.get("hgvs_column", defaults["hgvs_col"]),
            score_col=cfg.get("score_column", defaults["score_col"]),
            sd_col=cfg.get("sd_column", defaults["sd_col"]),
            output_path=snakemake.output[0],
        )
        return

    # Standalone CLI path
    parser = argparse.ArgumentParser(description="Format scores for MaveDB.")
    parser.add_argument("--scores", required=True, help="Score CSV from rosace or lilace")
    parser.add_argument("--output", required=True, help="Output MaveDB CSV path")
    parser.add_argument("--backend", default="rosace", choices=list(_BACKEND_DEFAULTS))
    parser.add_argument("--hgvs-column", default=None)
    parser.add_argument("--score-column", default=None)
    parser.add_argument("--sd-column", default=None)
    args = parser.parse_args()

    defaults = _BACKEND_DEFAULTS[args.backend]
    format_scores(
        scores_path=args.scores,
        backend=args.backend,
        hgvs_col=args.hgvs_column or defaults["hgvs_col"],
        score_col=args.score_column or defaults["score_col"],
        sd_col=args.sd_column or defaults["sd_col"],
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
