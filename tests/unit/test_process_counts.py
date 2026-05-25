# test_process_counts.py

import pytest
import os
import csv
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from workflow.rules.scripts.process_counts import (
    _run,
    process_sample,
    process_gatk_file,
)
from workflow.rules.scripts.script_utils import translate_orf


def _read_ref_aa(ref_dir, reference_fasta, orf_range="1-12"):
    """Helper: parse a reference FASTA and translate its ORF region.
    Mirrors what process_counts._run does once per rule invocation —
    tests that exercise process_sample directly need to do the same
    parsing themselves. Uses the shared translate_orf so the test path
    can't drift from production."""
    ref_path = os.path.join(ref_dir, reference_fasta)
    with open(ref_path) as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq
    return translate_orf(ref_sequence, orf_range)


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
    The output directory where process_sample writes files.
    """
    out_dir = tmp_path / "processed_counts"
    out_dir.mkdir()
    return str(out_dir)


@pytest.fixture
def mock_process_variants():
    """
    We'll patch the process_variants module to mock read_gatk_csv, process_variants_file,
    and write_enrich_df. That way we can isolate process_sample logic.
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


def test_process_sample_minimal(
    mock_ref_dir,
    reference_fasta,
    gatk_dir,
    output_dir,
    mock_process_variants,
):
    """
    Test that process_sample runs end-to-end for one sample, verifying
    expected files are created.
    """
    sample_name = "sampleA"
    gatk_csv = os.path.join(gatk_dir, f"{sample_name}.variantCounts")
    with open(gatk_csv, "w") as f:
        f.write("mock\ttsv\tlines")

    process_sample(
        sample_name=sample_name,
        experiment_name="experiment",
        ref_AA_sequence=_read_ref_aa(mock_ref_dir, os.path.basename(reference_fasta)),
        designed_df=pd.DataFrame(),
        max_deletion_length=3,
        noprocess=True,
        gatk_dir=gatk_dir,
        output_dir=output_dir,
    )

    # Check that process_sample wrote some expected files
    # 1) The processed CSV
    processed_csv = Path(output_dir) / f"{sample_name}.csv"
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
    lines = rejected_file.read_text().strip().split("\n")
    assert len(lines) <= 1, "Likely no rejected lines in the mock scenario"

    # 4) The stats files
    stats_dir = Path("stats") / "experiment" / "processing"
    assert stats_dir.exists(), "stats_dir should be created by process_sample"

    rejected_stats_file = stats_dir / f"{sample_name}_rejected_processing.tsv"
    accepted_stats_file = stats_dir / f"{sample_name}_accepted_processing.tsv"
    total_stats_file = stats_dir / f"{sample_name}_total_processing.tsv"

    assert rejected_stats_file.exists(), "Rejected stats should exist"
    assert accepted_stats_file.exists(), "Accepted stats should exist"
    assert total_stats_file.exists(), "Total stats should exist"


# -----------------------------------------------------------------------------
# _run-level tests for the per-sample rule entry point.
#
# After the per-sample split, _run handles exactly one sample identified by
# `snakemake.wildcards.sample_prefix`. The "missing variants_file under
# filtering" and "skip variants load under noprocess" guards still apply
# per-invocation; the cross-group caching test is no longer applicable
# (there's no loop to cache against).
# -----------------------------------------------------------------------------


def _build_run_snakemake(
    mock_snakemake, *, ref_dir, reference_fasta, variants_file,
    gatk_dir, noprocess, sample_prefix="sampleA",
):
    return mock_snakemake(
        config={
            "experiment": "experiment",
            "ref_dir": str(ref_dir),
            "reference": os.path.basename(reference_fasta),
            "orf": "1-12",
            "variants_file": str(variants_file),
            "noprocess": noprocess,
            "max_deletion_length": 3,
        },
        params={"gatk_dir": str(gatk_dir)},
        wildcards={"sample_prefix": sample_prefix},
        log=["/dev/null"],
    )


def test_run_raises_when_variants_file_missing_and_filtering(
    mock_snakemake, mock_ref_dir, reference_fasta, gatk_dir, output_dir,
):
    """Variant filtering needs a designed-variants CSV. If filtering is on
    (noprocess=False) but the file is missing, _run must fail loudly."""
    missing = Path(output_dir) / "does_not_exist.csv"
    snakemake = _build_run_snakemake(
        mock_snakemake,
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=missing,
        gatk_dir=gatk_dir,
        noprocess=False,
    )

    with pytest.raises(FileNotFoundError, match="Variants file not found"):
        _run(snakemake)


def test_run_skips_variants_file_load_when_noprocess(
    mock_snakemake, mock_ref_dir, reference_fasta, gatk_dir,
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
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=missing,
        gatk_dir=gatk_dir,
        noprocess=True,
        sample_prefix=sample_name,
    )

    # Should not raise.
    _run(snakemake)


def test_run_reads_inputs_once_per_invocation(
    mock_snakemake, mock_ref_dir, reference_fasta, gatk_dir, output_dir,
    mock_process_variants, mocker, tmp_path,
):
    """_run should parse the reference FASTA and read the variants CSV
    exactly once per invocation. (Each per-sample rule invocation is its
    own Python process, so this is a per-rule guarantee, not a cross-
    sample one.)"""
    # Real variants file so the read actually happens.
    variants_csv = tmp_path / "variants.csv"
    pd.DataFrame(
        [
            {"count": 0, "pos": 1, "mutation_type": "M", "name": "M1A",
             "codon": "GCT", "mutation": "A", "length": 1, "hgvs": "p.(M1A)"}
        ]
    ).to_csv(variants_csv, index=False)

    sample_name = "sampleA"
    with open(os.path.join(gatk_dir, f"{sample_name}.variantCounts"), "w") as f:
        f.write("mock\ttsv\tlines")

    snakemake = _build_run_snakemake(
        mock_snakemake,
        ref_dir=mock_ref_dir,
        reference_fasta=reference_fasta,
        variants_file=variants_csv,
        gatk_dir=gatk_dir,
        noprocess=False,
        sample_prefix=sample_name,
    )

    # Spy on the expensive operations.
    read_csv_spy = mocker.spy(pd, "read_csv")
    seqio_parse_spy = mocker.spy(SeqIO, "parse")

    _run(snakemake)

    # Reference is parsed exactly once.
    assert seqio_parse_spy.call_count == 1, (
        f"Expected SeqIO.parse called once; got {seqio_parse_spy.call_count}"
    )

    # pd.read_csv is called once for the variants CSV. (The experiment CSV
    # is no longer loaded per-sample; it's not needed when the wildcard
    # directly identifies the sample.)
    assert read_csv_spy.call_count == 1, (
        f"Expected pd.read_csv called once (variants only); "
        f"got {read_csv_spy.call_count}"
    )
