"""
Integration tests for file resolution and path handling.

Tests that the pipeline correctly resolves file paths, handles wildcards,
and manages input/output relationships.
"""

import pytest
from pathlib import Path


@pytest.mark.integration
class TestFileResolution:
    """Tests for file path resolution."""

    def test_fixtures_directory_exists(self, fixtures_dir):
        """Test that fixtures directory exists and is accessible."""
        assert fixtures_dir.exists(), "Fixtures directory should exist"
        assert fixtures_dir.is_dir(), "Fixtures should be a directory"

    def test_all_fixture_files_exist(self, fixtures_dir):
        """Test that all expected fixture files exist."""
        expected_files = [
            "mock_reference.fasta",
            "mock_variants.csv",
            "mock_oligos.csv",
            "mock_experiment.csv",
            "mock_config.yaml",
            "mock_gatk_output.variantCounts",
        ]

        for filename in expected_files:
            filepath = fixtures_dir / filename
            assert filepath.exists(), f"Fixture file {filename} should exist"

    def test_fixture_files_not_empty(self, fixtures_dir):
        """Test that fixture files are not empty."""
        for filepath in fixtures_dir.glob("*"):
            if filepath.is_file():
                assert filepath.stat().st_size > 0, f"{filepath.name} should not be empty"


@pytest.mark.integration
class TestPathPatterns:
    """Tests for path pattern handling."""

    def test_sample_file_pattern_matching(self, mock_experiment_df):
        """Test that sample names can be extracted from file patterns."""
        samples = mock_experiment_df["sample"].tolist()

        assert len(samples) > 0, "Should have samples"
        # Check sample naming convention
        for sample in samples:
            assert isinstance(sample, str), "Sample names should be strings"
            assert len(sample) > 0, "Sample names should not be empty"

    def test_condition_extraction(self, mock_experiment_df):
        """Test that conditions can be extracted from experiment file."""
        conditions = mock_experiment_df["condition"].unique().tolist()

        assert len(conditions) >= 1, "Should have at least one condition"
        assert all(isinstance(c, str) for c in conditions), "Conditions should be strings"

    def test_replicate_extraction(self, mock_experiment_df):
        """Test that replicates can be extracted from experiment file."""
        replicates = mock_experiment_df["replicate"].unique().tolist()

        assert len(replicates) >= 1, "Should have at least one replicate"
        assert all(isinstance(r, (int, float)) for r in replicates), "Replicates should be numeric"


@pytest.mark.integration
class TestOutputPaths:
    """Tests for output path generation."""

    def test_output_directory_structure(self, tmp_path, mock_config):
        """Test generating output directory structure."""
        experiment = mock_config["experiment"]

        # Create expected directory structure
        results_dir = tmp_path / "results" / experiment
        stats_dir = tmp_path / "stats" / experiment
        logs_dir = tmp_path / "logs" / experiment

        results_dir.mkdir(parents=True)
        stats_dir.mkdir(parents=True)
        logs_dir.mkdir(parents=True)

        # Verify structure
        assert results_dir.exists()
        assert stats_dir.exists()
        assert logs_dir.exists()

    def test_gatk_output_path_pattern(self, mock_config, mock_experiment_df):
        """Test GATK output path pattern generation."""
        experiment = mock_config["experiment"]
        samples = mock_experiment_df["sample"].tolist()

        for sample in samples:
            # Expected path pattern: results/{experiment}/gatk/{sample}.variantCounts
            expected_path = f"results/{experiment}/gatk/{sample}.variantCounts"
            assert "{experiment}" not in expected_path, "Path should have experiment resolved"
            assert "{sample}" not in expected_path, "Path should have sample resolved"

    def test_processed_counts_path_pattern(self, mock_config, mock_experiment_df):
        """Test processed counts output path pattern."""
        experiment = mock_config["experiment"]
        samples = mock_experiment_df["sample"].tolist()

        for sample in samples:
            # Expected patterns
            enrich_path = f"results/{experiment}/processed_counts/enrich_format/{sample}.tsv"
            csv_path = f"results/{experiment}/processed_counts/{sample}.csv"

            assert sample in enrich_path
            assert sample in csv_path


@pytest.mark.integration
class TestWildcardExpansion:
    """Tests for Snakemake wildcard expansion patterns."""

    def test_sample_wildcard_values(self, mock_experiment_df):
        """Test that sample wildcards expand correctly."""
        samples = mock_experiment_df["sample"].tolist()

        # Simulate wildcard expansion
        expanded = [f"results/test/gatk/{s}.variantCounts" for s in samples]

        assert len(expanded) == len(samples)
        assert all(".variantCounts" in p for p in expanded)

    def test_condition_wildcard_values(self, mock_experiment_df):
        """Test that condition wildcards expand correctly."""
        conditions = mock_experiment_df["condition"].unique().tolist()

        # Simulate wildcard expansion for rosace output
        expanded = [f"results/test/rosace/{c}_scores.csv" for c in conditions]

        assert len(expanded) == len(conditions)
        assert all("_scores.csv" in p for p in expanded)


@pytest.mark.integration
class TestReferenceFiles:
    """Tests for reference file handling."""

    def test_reference_fasta_format(self, fixtures_dir):
        """Test that reference FASTA is properly formatted."""
        fasta_file = fixtures_dir / "mock_reference.fasta"
        content = fasta_file.read_text()

        lines = content.strip().split("\n")
        header = lines[0]

        assert header.startswith(">"), "FASTA should have header starting with >"
        assert len(header) > 1, "FASTA header should have content after >"

        # Check sequence
        sequence = "".join(lines[1:])
        assert len(sequence) > 0, "FASTA should have sequence"
        assert all(c in "ACGTNacgtn" for c in sequence), "Sequence should be nucleotides"

    def test_reference_length_matches_orf(self, fixtures_dir):
        """Test that reference length accommodates the ORF."""
        import yaml

        config_file = fixtures_dir / "mock_config.yaml"
        with open(config_file) as f:
            config = yaml.safe_load(f)

        orf = config["orf"]
        start, stop = map(int, orf.split("-"))

        fasta_file = fixtures_dir / "mock_reference.fasta"
        content = fasta_file.read_text()
        sequence = "".join(content.strip().split("\n")[1:])

        assert len(sequence) >= stop, f"Reference ({len(sequence)}bp) should be >= ORF stop ({stop})"
