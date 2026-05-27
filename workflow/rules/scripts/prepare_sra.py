"""Prepare SRA submission metadata and file list.

Generates:
  - sra_metadata.tsv   SRA run metadata spreadsheet (ready for NCBI upload)
  - sra_files.txt      Paths of all FASTQ files to upload

Run via Snakemake or standalone:
  python prepare_sra.py --experiment-file <exp.csv> --data-dir <dir/> --output-dir <out/>
"""

import argparse
import re
from pathlib import Path

import pandas as pd

_R1_PAT = re.compile(r"[._](?:R1|1)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")
_R2_PAT = re.compile(r"[._](?:R2|2)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")

# SRA requires a controlled vocabulary for instrument model.
# User must fill this in; listed here as reminder.
_INSTRUMENT_PLACEHOLDER = "FILL_IN: e.g. Illumina NovaSeq 6000"
_BIOSAMPLE_PLACEHOLDER = "FILL_IN: e.g. SAMN00000000"


def resolve_fastq_pair(data_dir: Path, prefix: str):
    R1, R2 = None, None
    for f in data_dir.glob(f"{prefix}*"):
        if _R1_PAT.search(f.name) and R1 is None:
            R1 = f
        elif _R2_PAT.search(f.name) and R2 is None:
            R2 = f
    return R1, R2


def build_metadata(experiments: pd.DataFrame, data_dir: Path, experiment_name: str, sra_cfg: dict):
    instrument = sra_cfg.get("instrument_model", _INSTRUMENT_PLACEHOLDER)
    design_desc = sra_cfg.get(
        "design_description",
        f"Deep mutational scanning amplicon sequencing for experiment {experiment_name}",
    )

    rows = []
    all_files = []
    seen_prefixes = {}  # prefix -> row index, for deduplication warnings

    for _, row in experiments.iterrows():
        sample = row["sample"]
        condition = row["condition"]
        replicate = row["replicate"]
        prefix = row["file"]

        R1, R2 = resolve_fastq_pair(data_dir, prefix)

        r1_name = R1.name if R1 else f"NOT_FOUND ({prefix} R1)"
        r2_name = R2.name if R2 else f"NOT_FOUND ({prefix} R2)"

        if prefix in seen_prefixes:
            # Multiple samples point to the same FASTQ — note both in the title
            prev_idx = seen_prefixes[prefix]
            prev_title = rows[prev_idx]["title"]
            rows[prev_idx]["title"] = prev_title + f"; also {condition} replicate {replicate}"
            if R1:
                all_files.append(str(R1))
            if R2:
                all_files.append(str(R2))
            continue

        seen_prefixes[prefix] = len(rows)
        rows.append(
            {
                "sample_name": sample,
                "library_ID": prefix,
                "title": f"{experiment_name}: {condition} replicate {replicate}",
                "library_strategy": "AMPLICON",
                "library_source": "GENOMIC",
                "library_selection": "PCR",
                "library_layout": "paired",
                "platform": "ILLUMINA",
                "instrument_model": instrument,
                "design_description": design_desc,
                "filetype": "fastq",
                "filename": r1_name,
                "filename2": r2_name,
                "biosample_accession": _BIOSAMPLE_PLACEHOLDER,
            }
        )

        if R1:
            all_files.append(str(R1))
        if R2:
            all_files.append(str(R2))

    return pd.DataFrame(rows), sorted(set(all_files))


def main():
    # Snakemake execution path. `"snakemake" in dir()` would always
    # return False inside a function — dir() returns local names, not
    # module globals — so use globals() to detect the snakemake-injected
    # variable. Same fix as format_mavedb.py (2026-05-27).
    if "snakemake" in globals():
        experiments = pd.read_csv(snakemake.input.experiment_file)
        data_dir = Path(snakemake.params.data_dir).resolve()
        sra_cfg = snakemake.params.sra_config or {}
        experiment_name = snakemake.config["experiment"]

        metadata, all_files = build_metadata(experiments, data_dir, experiment_name, sra_cfg)
        metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)
        with open(snakemake.output.filelist, "w") as fh:
            fh.write("\n".join(all_files) + ("\n" if all_files else ""))
        return

    # Standalone CLI path
    parser = argparse.ArgumentParser(description="Prepare SRA submission package.")
    parser.add_argument("--experiment-file", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--experiment-name", default="experiment")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--instrument-model", default=None)
    parser.add_argument("--design-description", default=None)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sra_cfg = {}
    if args.instrument_model:
        sra_cfg["instrument_model"] = args.instrument_model
    if args.design_description:
        sra_cfg["design_description"] = args.design_description

    experiments = pd.read_csv(args.experiment_file)
    data_dir = Path(args.data_dir).resolve()

    metadata, all_files = build_metadata(experiments, data_dir, args.experiment_name, sra_cfg)
    metadata.to_csv(out_dir / "sra_metadata.tsv", sep="\t", index=False)
    with open(out_dir / "sra_files.txt", "w") as fh:
        fh.write("\n".join(all_files) + ("\n" if all_files else ""))

    print(f"Wrote {len(metadata)} SRA run rows to {out_dir}/sra_metadata.tsv")
    print(f"Listed {len(all_files)} FASTQ files in {out_dir}/sra_files.txt")


if __name__ == "__main__":
    main()
