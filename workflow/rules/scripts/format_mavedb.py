"""Format variant scores into MaveDB score set CSV format.

Output columns: hgvs_pro, score, sd (if available)
Run via Snakemake or standalone:
  python format_mavedb.py --scores <scores.csv> --output <out.csv> [--backend rosace]
"""

import argparse
import re
import sys

import pandas as pd

# Default column names per scoring backend.
# Override via config: mavedb.hgvs_column, mavedb.score_column, mavedb.sd_column
_BACKEND_DEFAULTS = {
    "rosace": {"hgvs_col": "variants", "score_col": "mean", "sd_col": "sd"},
    "lilace": {"hgvs_col": "variants", "score_col": "mean", "sd_col": "sd"},
}

_SYNONYMOUS_RE = re.compile(r"^p\.\(([A-Z])(\d+)([A-Z])\)$")


def normalize_hgvs(hgvs: str) -> str:
    """Convert same-AA synonymous notation p.(X1X) → p.(X1=) per HGVS standard."""
    m = _SYNONYMOUS_RE.match(hgvs)
    if m and m.group(1) == m.group(3):
        return f"p.({m.group(1)}{m.group(2)}=)"
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
    out["hgvs_pro"] = scores[hgvs_col].apply(normalize_hgvs)
    out["score"] = scores[score_col]

    if sd_col and sd_col in scores.columns:
        out["sd"] = scores[sd_col]

    out.to_csv(output_path, index=False)


def main():
    # Snakemake execution path
    if "snakemake" in dir():
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


main()
