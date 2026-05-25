# test_process_counts.py

import pytest
import os
import csv
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from workflow.rules.scripts.process_counts import (
    _run,
    process_experiment,
    process_gatk_file,
)


def _read_ref_aa(ref_dir, reference_fasta, orf_range="1-12"):
    """Helper: parse a reference FASTA and translate its ORF region.
    Mirrors what process_counts._run does once per pipeline invocation —
    tests that exercise process_experiment directly need to do the same
    parsing themselves now that the responsibility moved out of the
    per-call function (audit item M7)."""
    ref_path = os.path.join(ref_dir, reference_fasta)
    with open(ref_path) as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq
    start, end = map(int, orf_range.split("-"))
    return ref_sequence[start - 1 : end].translate()


@pytest.fixture
def mock_ref_dir(tmp_path):
    """
    Create a temporary directory that will contain a mock reference FASTA file.
    """
    ref_dir = tmp_path / "ref_dir"
    ref_dir.mkdir()
    return ref_dir


@pytest.fixture
def reference_fasta(mock_ref_dir):
    """
    Create a small reference FASTA in the mock ref_dir, returning the filename.
    We'll create a single record with a short sequence so we can test orf extraction.
    """
    fasta_file = mock_ref_dir / "test_ref.fasta"
    record = SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), id="test_ref")
    SeqIO.write(record, fasta_file, "fasta")
    return str(fasta_file)


@pytest.fixture
def variants_file(tmp_path):
    """
    Create a small 'designed variants' CSV file with the columns:
      count, pos, mutation_type, name, codon, mutation, length, hgvs
    We'll place it in the tmp_path and return its path.
    """
    variants_csv = tmp_path / "designed_variants.csv"
    data = [
        # Header
        [
            "count",
            "pos",
            "mutation_type",
            "name",
            "codon",
            "mutation",
            "length",
            "hgvs",
        ],
        # Rows
        [0, 1, "M", "M1A", "AAA", "A", 1, "p.(M1A)"],
        [0, 2, "S", "S2G", "", "G", 1, "p.(S2G)"],
    ]
    with variants_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(data)

    return str(variants_csv)


@pytest.fixture
def oligo_file(tmp_path):
    """
    If you're testing the 'regenerate_variants' logic, create a mock oligo file
    for generate_variants. If not, we can skip this or provide an empty file.
    """
    oligo_file = tmp_path / "oligos.txt"
    oligo_file.write_text("Mock oligo data\n")
    return str(oligo_file)


@pytest.fixture
def gatk_dir(tmp_path):
    """
    A directory to store mock GATK CSV outputs (one per sample).
    We'll write minimal data for each sample.
    """
    gatk_path = tmp_path / "gatk"
    gatk_path.mkdir()
    return str(gatk_path)


@pytest.fixture
def output_dir(tmp_path):
    """
    The output directory where process_experiment writes files.
    """
    out_dir = tmp_path / "processed_counts"
    out_dir.mkdir()
    return str(out_dir)


@pytest.fixture
def mock_process_variants():
    """
    We'll patch the process_variants module to mock read_gatk_csv, process_variants_file,
    and write_enrich_df. That way we can isolate process_experiment logic.
    """
    with patch("workflow.rules.scripts.process_counts.process_variants") as mock_mod:
        # read_gatk_csv will return a small list of lines
        mock_mod.read_gatk_csv.return_value = [
            ["5", "0", "0", "3", "A", "1", "1:AAA>CCC", "Mstuff", "M1A"],
            ["2", "0", "0", "3", "A", "1", "2:AAA>AAA", "Sstuff", ""],
        ]

        # process_variants_file returns df, other, rejected_stats, accepted_stats, total_stats
        # We'll return a mock DataFrame, an empty 'other', and some stats
        def mock_process_variants_file(
            gatk_list, designed_df, ref_AA_seq, max_del_len, noprocess
        ):
            df = pd.DataFrame(
                {
                    "hgvs": ["p.(M1A)", "p.(S2G)"],
                    "count": [5, 2],
                    "mutation_type": ["M", "S"],
                }
            )
            other = []
            rejected_stats = {"fs_counts": 0}
            accepted_stats = {"accepted_sub_counts": 5, "accepted_syn_counts": 2}
            total_stats = {"total_counts": 7}
            return df, other, rejected_stats, accepted_stats, total_stats

        mock_mod.process_variants_file.side_effect = mock_process_variants_file

        # write_enrich_df writes a minimal enrich2-like file for assertions.
        def mock_write_enrich_df(file_path, _df, _noprocess):
            p = Path(file_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("hgvs\tcount\np.(M1A)\t5\n")

        mock_mod.write_enrich_df.side_effect = mock_write_enrich_df

        # write_stats_file writes minimal content for assertions.
        def mock_write_stats_file(file_path, _stats):
            p = Path(file_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("Mutation\tCount\n")

        mock_mod.write_stats_file.side_effect = mock_write_stats_file

        yield mock_mod


def test_process_experiment_minimal(
    mock_ref_dir,
    reference_fasta,
    gatk_dir,
    output_dir,
    mock_process_variants,
):
    """
    Test that process_experiment runs end-to-end for a minimal set of samples,
    verifying certain files are created and we see expected logs.
    """
    # We'll create a minimal GATK CSV for one sample named 'sampleA'
    sample_name = "sampleA"
    gatk_csv = os.path.join(gatk_dir, f"{sample_name}.variantCounts")
    with open(gatk_csv, "w") as f:
        f.write("mock\ttsv\tlines")

    process_experiment(
        sample_list=[sample_name],
        experiment_name="experiment",
        ref_AA_sequence=_read_ref_aa(mock_ref_dir, os.path.basename(reference_fasta)),
        designed_df=pd.DataFrame(),
        max_deletion_length=3,
        noprocess=True,
        gatk_dir=gatk_dir,
        output_dir=output_dir,
    )

    # Check that process_experiment wrote some expected files
    # 1) The processed CSV
    processed_csv = Path(output_dir) / f"{sample_name}.csv"
    # Print the contents of dir
    print(list(Path(output_dir).rglob("*")))
    assert processed_csv.exists(), "Processed CSV should be created"
    df_processed = pd.read_csv(processed_csv)
    # We expect the patched process_variants_file returned 2 rows:
    assert len(df_processed) == 2

    # 2) The enrich_format file
    enrich_file = (
        Path(output_dir) / "enrich_format" / f"{sample_name}.tsv"
    )
    assert enrich_file.exists(), "Enrich2-readable file should be created"

    # 3) The rejected variants
    rejected_file = Path(output_dir) / "rejected" / f"rejected_{sample_name}.csv"
    assert rejected_file.exists(), "Rejected variants file should be created"
    # Possibly empty, since we mocked a minimal scenario
    lines = rejected_file.read_text().strip().split("\n")
    assert len(lines) <= 1, "Likely no rejected lines in the mock scenario"

    # 4) The stats files
    stats_dir = Path("stats") / "experiment" / "processing"
    assert stats_dir.exists(), "stats_dir should be created by process_experiment"

    # We'll check for sampleA_rejected_processing.tsv, etc
    rejected_stats_file = stats_dir / f"{sample_name}_rejected_processing.tsv"
    accepted_stats_file = stats_dir / f"{sample_name}_accepted_processing.tsv"
    total_stats_file = stats_dir / f"{sample_name}_total_processing.tsv"

    assert rejected_stats_file.exists(), "Rejected stats should exist"
    assert accepted_stats_file.exists(), "Accepted stats should exist"
    assert total_stats_file.exists(), "Total stats should exist"
    # Quick check of contents
    content = rejected_stats_file.read_text()
    # We'll not do detailed checks here, as we used a patched version.
    # But you can parse or do line-by-line checks if needed.


# -----------------------------------------------------------------------------
# _run-level tests for behavior hoisted out of process_experiment (M7).
#
# The "missing variants_file under filtering" check and the "skip variants
# load under noprocess" shortcut both used to live inside process_experiment.
# They moved into _run when the reference + variants-file reads were hoisted
# out of the per-replicate-group loop. The tests follow.
# -----------------------------------------------------------------------------


@pytest.fixture
def experiment_csv(tmp_path):
    """Tiny experiment CSV with one sample, one condition, one replicate."""
    csv_file = tmp_path / "experiment.csv"
    pd.DataFrame(
        [{"sample": "sampleA", "condition": "A", "replicate": 1, "time": 0}]
    ).to_csv(csv_file, index=False)
    return csv_file


def _build_run_snakemake(
    mock_snakemake, *, experiment_csv, ref_dir, reference_fasta, variants_file,
    gatk_dir, noprocess, tiled=False,
):
    return mock_snakemake(
        config={
            "experiment": "experiment",
            "experiment_file": str(experiment_csv),
            "ref_dir": str(ref_dir),
            "reference": os.path.basename(reference_fasta),
            "orf": "1-12",
            "variants_file": str(variants_file),
            "noprocess": noprocess,
            "max_deletion_length": 3,
            "tiled": tiled,
        },
        params={"gatk_dir": str(gatk_dir)},
        log=["/dev/null"],
    )


def test_run_raises_when_variants_file_missing_and_filtering(
    mock_snakemake, mock_ref_dir, reference_fasta, experiment_csv, gatk_dir, output_dir,
):
    """Variant filtering needs a designed-variants CSV. If filtering is on
    (noprocess=False) but the file is missing, _run must fail loudly."""
    missing = Path(output_dir) / "does_not_exist.csv"
    snakemake = _build_run_snakemake(
        mock_snakemake,
        experiment_csv=experiment_csv,
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=missing,
        gatk_dir=gatk_dir,
        noprocess=False,
    )

    with pytest.raises(FileNotFoundError, match="Variants file not found"):
        _run(snakemake)


def test_run_skips_variants_file_load_when_noprocess(
    mock_snakemake, mock_ref_dir, reference_fasta, experiment_csv, gatk_dir,
    output_dir, mock_process_variants,
):
    """noprocess=True must NOT require the variants CSV to exist — the
    file path can be invalid and _run should still complete."""
    missing = Path(output_dir) / "does_not_exist.csv"

    # Write a minimal GATK file that the (mocked) process_variants will read.
    sample_name = "sampleA"
    gatk_csv = os.path.join(gatk_dir, f"{sample_name}.variantCounts")
    with open(gatk_csv, "w") as f:
        f.write("mock\ttsv\tlines")

    snakemake = _build_run_snakemake(
        mock_snakemake,
        experiment_csv=experiment_csv,
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=missing,
        gatk_dir=gatk_dir,
        noprocess=True,
    )

    # Should not raise.
    _run(snakemake)

    # The per-sample outputs that don't depend on the variants file should
    # still have been produced. Note: _run writes outputs under
    # results/{experiment_name}/processed_counts/ in the current working
    # directory, not under our test output_dir fixture. Don't assert
    # specific paths — just that _run returned without raising.


def test_run_loads_reference_and_variants_once_regardless_of_groups(
    mock_snakemake, mock_ref_dir, reference_fasta, gatk_dir, output_dir,
    mock_process_variants, mocker, tmp_path,
):
    """M7 regression guard. The reference FASTA parse + ORF translate and
    the variants-file read should happen once per _run invocation, not
    once per (condition × tile × replicate) call to process_experiment."""
    # Build an experiment CSV with multiple replicate groups so the old
    # implementation would re-read everything per group.
    csv_file = tmp_path / "multi_group_experiment.csv"
    pd.DataFrame(
        [
            {"sample": "A_T0", "condition": "A", "replicate": 1, "time": 0},
            {"sample": "A_T1", "condition": "A", "replicate": 1, "time": 1},
            {"sample": "A_T0_r2", "condition": "A", "replicate": 2, "time": 0},
            {"sample": "A_T1_r2", "condition": "A", "replicate": 2, "time": 1},
            {"sample": "B_T0", "condition": "B", "replicate": 1, "time": 0},
            {"sample": "B_T1", "condition": "B", "replicate": 1, "time": 1},
        ]
    ).to_csv(csv_file, index=False)

    # Real variants file (small) so the read actually happens.
    variants_csv = tmp_path / "variants.csv"
    pd.DataFrame(
        [
            {"count": 0, "pos": 1, "mutation_type": "M", "name": "M1A",
             "codon": "GCT", "mutation": "A", "length": 1, "hgvs": "p.(M1A)"}
        ]
    ).to_csv(variants_csv, index=False)

    # Write a GATK file per sample so each replicate group has something
    # to process.
    for name in ("A_T0", "A_T1", "A_T0_r2", "A_T1_r2", "B_T0", "B_T1"):
        with open(os.path.join(gatk_dir, f"{name}.variantCounts"), "w") as f:
            f.write("mock\ttsv\tlines")

    snakemake = _build_run_snakemake(
        mock_snakemake,
        experiment_csv=csv_file,
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=variants_csv,
        gatk_dir=gatk_dir,
        noprocess=False,
    )

    # Spy on the expensive operations.
    read_csv_spy = mocker.spy(pd, "read_csv")
    seqio_parse_spy = mocker.spy(SeqIO, "parse")

    _run(snakemake)

    # Reference is parsed exactly once.
    assert seqio_parse_spy.call_count == 1, (
        f"Expected SeqIO.parse called once; got {seqio_parse_spy.call_count}"
    )

    # pd.read_csv is called once for the experiment CSV and once for the
    # variants CSV — total 2, regardless of replicate-group count. The
    # buggy pre-M7 version would call it (1 + N) times where N is the
    # number of replicate groups (3 here: A/1, A/2, B/1).
    assert read_csv_spy.call_count == 2, (
        f"Expected pd.read_csv called twice (experiment + variants); "
        f"got {read_csv_spy.call_count}"
    )
