"""
Tests for generate_baseline_file_list._run.

The historical bugs (audit item M1):

1. FastQC paths were built by concatenating the *sample name* with `R1_001`
   etc., but FastQC names its outputs from the actual *fastq filename*. The
   pipeline already computes that mapping (via common.smk's `fastqc_names`
   dict, fed by `resolve_fastq_pair`); the script needs to consult it
   rather than guessing.

2. The processed-counts path was `processed_counts/counts/{sample}.csv`,
   but `process_counts.py` writes to `processed_counts/{sample}.csv`
   (no `counts/` subdirectory).

Both meant MultiQC silently skipped those entries from the baseline
report. Tests cover the corrected paths.
"""

import pandas as pd

from workflow.rules.scripts.generate_baseline_file_list import _run


def _write_experiment_csv(tmp_path):
    """Two baseline samples sharing one fastq file prefix `mock_reads`."""
    rows = [
        {"sample": "B_T0", "condition": "baseline", "replicate": 1, "time": 0,
         "tile": 1, "file": "mock_reads"},
        {"sample": "B_T1", "condition": "baseline", "replicate": 1, "time": 1,
         "tile": 1, "file": "mock_reads"},
        # An experimental (non-baseline) sample that must NOT appear in the
        # baseline file list.
        {"sample": "A_T0", "condition": "treated", "replicate": 1, "time": 0,
         "tile": 1, "file": "mock_reads"},
    ]
    csv = tmp_path / "experiment.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    return csv


def test_run_writes_fastqc_paths_with_resolved_names(tmp_path, mock_snakemake):
    """FastQC paths must use the names from the fastqc_names map, not the
    sample name from the experiment CSV."""
    experiment_csv = _write_experiment_csv(tmp_path)
    output_file = tmp_path / "baseline_file_list.txt"

    fastqc_names = {
        "mock_reads": {
            "R1": "mock_reads_R1_001",
            "R2": "mock_reads_R2_001",
        },
    }

    snakemake = mock_snakemake(
        config={
            "experiment": "expt",
            "experiment_file": str(experiment_csv),
            "baseline_condition": "baseline",
        },
        params={"fastqc_names": fastqc_names},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    lines = output_file.read_text().splitlines()

    # FastQC: must be named off the fastq file's prefix, not the sample name.
    assert "./stats/expt/fastqc/mock_reads_R1_001_fastqc.html" in lines
    assert "./stats/expt/fastqc/mock_reads_R2_001_fastqc.html" in lines

    # And NOT the old buggy shape (sample name + "R1_001", no separator).
    assert not any("B_T0R1_001_fastqc.html" in line for line in lines)
    assert not any("B_T1R1_001_fastqc.html" in line for line in lines)


def test_run_writes_processed_counts_without_counts_subdir(tmp_path, mock_snakemake):
    """processed_counts CSVs live directly under processed_counts/, not
    under a `counts/` subdirectory that doesn't exist."""
    experiment_csv = _write_experiment_csv(tmp_path)
    output_file = tmp_path / "baseline_file_list.txt"

    snakemake = mock_snakemake(
        config={
            "experiment": "expt",
            "experiment_file": str(experiment_csv),
            "baseline_condition": "baseline",
        },
        params={"fastqc_names": {"mock_reads": {"R1": "x", "R2": "y"}}},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    text = output_file.read_text()

    assert "./results/expt/processed_counts/B_T0.csv" in text
    assert "./results/expt/processed_counts/B_T1.csv" in text
    assert "/processed_counts/counts/" not in text, (
        "processed_counts/counts/ path should be gone — process_counts.py "
        "writes directly to processed_counts/"
    )


def test_run_excludes_non_baseline_samples(tmp_path, mock_snakemake):
    """Only samples with condition == baseline_condition should appear in the
    file list. Confirms the rule's intent isn't accidentally pulling in
    treated samples."""
    experiment_csv = _write_experiment_csv(tmp_path)
    output_file = tmp_path / "baseline_file_list.txt"

    snakemake = mock_snakemake(
        config={
            "experiment": "expt",
            "experiment_file": str(experiment_csv),
            "baseline_condition": "baseline",
        },
        params={"fastqc_names": {"mock_reads": {"R1": "r1", "R2": "r2"}}},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    text = output_file.read_text()
    # B_T0 and B_T1 are baseline -> present.
    assert "B_T0" in text
    assert "B_T1" in text
    # A_T0 is treated -> absent.
    assert "A_T0" not in text
