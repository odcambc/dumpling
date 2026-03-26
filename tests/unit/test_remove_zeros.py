import pytest
import pandas as pd
import os
from pathlib import Path

# Import the functions you want to test directly
# Adjust if the script is named differently or in a different location
from workflow.rules.scripts.remove_zeros import remove_zeros_enrich, remove_zeros


@pytest.fixture
def mock_data_dir(tmp_path):
    """
    Create a temporary directory structure with mock TSV files for testing.
    We'll place a few 'enrich_format' TSV files that have a 'hgvs' column
    and a 'count' column for each sample.
    """
    data_dir = tmp_path / "mock_data"
    data_dir.mkdir()

    # Example contents for two sample files
    # File 1
    sample1_file = data_dir / "sample1.tsv"
    sample1_content = "hgvs\tcount\n" "varA\t10\n" "varB\t0\n" "varC\t5\n"
    sample1_file.write_text(sample1_content)

    # File 2
    sample2_file = data_dir / "sample2.tsv"
    sample2_content = "hgvs\tcount\n" "varA\t0\n" "varB\t0\n" "varC\t5\n" "varD\t0\n"
    sample2_file.write_text(sample2_content)

    # Return the directory path
    return data_dir


def test_remove_zeros_enrich(tmp_path, mock_data_dir):
    """
    Test remove_zeros_enrich with two input files:
      - sample1 has varA=10, varB=0, varC=5
      - sample2 has varA=0, varB=0, varC=5, varD=0

    'varB' is 0 in both, 'varD' is 0 in sample2 and absent in sample1 => also 0 for sample1
    => 'varB' and 'varD' are unobserved across *all* files, so they should be removed.
    """
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    sample1 = (mock_data_dir / "sample1.tsv").as_posix()
    sample2 = (mock_data_dir / "sample2.tsv").as_posix()
    file_list = [sample1, sample2]

    unobserved = remove_zeros_enrich(file_list, output_dir.as_posix())

    # We expect unobserved variants => varB, varD
    assert set(unobserved) == {"varB", "varD"}

    # Check the rewritten files in output_dir
    # sample1.tsv
    out1 = output_dir / "sample1.tsv"
    assert out1.exists()
    df1 = pd.read_csv(out1, sep="\t")
    # This should only have varA, varC
    # varA=10, varC=5
    assert len(df1) == 2
    assert set(df1["hgvs"].tolist()) == {"varA", "varC"}
    # sample2.tsv
    out2 = output_dir / "sample2.tsv"
    df2 = pd.read_csv(out2, sep="\t")
    assert len(df2) == 2
    assert set(df2["hgvs"].tolist()) == {"varA", "varC"}

    # varA in sample2 is 0 originally, but not removed because sample1 had varA=10
    # Here in the filtered output, varA is present but should be 0 in sample2:
    row_varA = df2.loc[df2["hgvs"] == "varA", "count"].iloc[0]
    assert row_varA == 0  # because sample2 had varA=0
    row_varC = df2.loc[df2["hgvs"] == "varC", "count"].iloc[0]
    assert row_varC == 5


def test_remove_zeros(tmp_path, mock_data_dir):
    """
    Test the remove_zeros function, which calls remove_zeros_enrich internally
    and writes an 'unobserved_variants.csv' file.

    We'll pass a list of sample names: e.g., ['sample1', 'sample2'].
    We'll check that the 'rejected' folder has the unobserved file
    with the correct variant set.
    """
    output_dir = tmp_path / "filtered"
    output_dir.mkdir(parents=True)

    # We pass an "input_dir" that matches mock_data_dir
    experiment_list = ["sample1", "sample2"]
    input_dir = f"{mock_data_dir.as_posix()}/"
    experiment_name = "mock_experiment"

    remove_zeros(experiment_list, input_dir, output_dir.as_posix(), experiment_name)

    # Check output files
    filtered1 = output_dir / "sample1.tsv"
    assert filtered1.exists()
    filtered2 = output_dir / "sample2.tsv"
    assert filtered2.exists()

    # Check the 'rejected' folder for unobserved variants
    rejected_dir = output_dir / "rejected"
    assert rejected_dir.exists()

    unobserved_file = rejected_dir / "mock_experiment_unobserved_variants.csv"
    assert unobserved_file.exists()
    lines = unobserved_file.read_text().strip().split("\n")
    # We expect "varB" and "varD" as in previous test
    assert set(lines) == {"varB", "varD"}
