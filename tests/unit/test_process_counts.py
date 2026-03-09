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
    process_experiment,
    write_counts_output,
    write_stats_output,
    process_gatk_file,
)


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
    for process_oligo_list. If not, we can skip this or provide an empty file.
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
    oligo_file,
    variants_file,
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
        ref_dir=mock_ref_dir,
        reference_fasta=os.path.basename(
            reference_fasta
        ),  # since we pass ref_dir separately
        oligo_file=oligo_file,
        variants_file=variants_file,
        orf_range="1-12",  # We'll just parse 12 bases => 4 amino acids
        max_deletion_length=3,
        noprocess=True,
        gatk_dir=gatk_dir,
        output_dir=output_dir,
        regenerate_variants=True,
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
