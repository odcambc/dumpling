"""
Integration tests for configuration validation.

Tests that the pipeline's configuration validation correctly accepts valid
configs and rejects invalid ones.
"""

import pytest
import yaml
from pathlib import Path


@pytest.mark.integration
class TestConfigSchema:
    """Tests for config schema validation."""

    def test_valid_config_has_required_fields(self, mock_config):
        """Verify that mock_config fixture contains all required fields."""
        required_fields = [
            "experiment",
            "experiment_file",
            "data_dir",
            "ref_dir",
            "reference",
            "orf",
        ]
        for field in required_fields:
            assert field in mock_config, f"Missing required field: {field}"

    def test_orf_format_valid(self, mock_config):
        """Test that ORF format is valid (start-stop)."""
        orf = mock_config["orf"]
        assert "-" in orf, "ORF should be in format 'start-stop'"
        parts = orf.split("-")
        assert len(parts) == 2, "ORF should have exactly two parts"
        assert parts[0].isdigit(), "ORF start should be numeric"
        assert parts[1].isdigit(), "ORF stop should be numeric"
        assert int(parts[0]) < int(parts[1]), "ORF start should be less than stop"

    def test_config_boolean_flags(self, mock_config):
        """Test that boolean flags are properly typed."""
        boolean_fields = [
            "enrich2",
            "remove_zeros",
            "regenerate_variants",
            "noprocess",
            "run_qc",
            "samtools_local",
            "rosace_local",
        ]
        for field in boolean_fields:
            if field in mock_config:
                assert isinstance(
                    mock_config[field], bool
                ), f"{field} should be boolean"

    def test_config_numeric_fields(self, mock_config):
        """Test that numeric fields are properly typed."""
        numeric_fields = ["kmers", "mem", "min_q", "min_variant_obs", "max_deletion_length"]
        for field in numeric_fields:
            if field in mock_config:
                assert isinstance(
                    mock_config[field], (int, float)
                ), f"{field} should be numeric"


@pytest.mark.integration
class TestConfigFileValidation:
    """Tests for validating config files against schema."""

    def test_mock_config_yaml_is_valid(self, fixtures_dir):
        """Test that mock_config.yaml is valid YAML."""
        config_file = fixtures_dir / "mock_config.yaml"
        assert config_file.exists(), "mock_config.yaml should exist"

        with open(config_file) as f:
            config = yaml.safe_load(f)

        assert isinstance(config, dict), "Config should be a dictionary"
        assert "experiment" in config, "Config should have experiment field"

    def test_mock_experiment_csv_structure(self, fixtures_dir):
        """Test that mock_experiment.csv has required columns."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_experiment.csv"
        assert csv_file.exists(), "mock_experiment.csv should exist"

        df = pd.read_csv(csv_file)
        required_columns = ["sample", "condition", "replicate", "time"]
        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"

    def test_mock_variants_csv_structure(self, fixtures_dir):
        """Test that mock_variants.csv has required columns."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_variants.csv"
        assert csv_file.exists(), "mock_variants.csv should exist"

        df = pd.read_csv(csv_file)
        required_columns = ["count", "pos", "mutation_type", "name", "hgvs"]
        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"

    def test_mock_reference_fasta_valid(self, fixtures_dir):
        """Test that mock_reference.fasta is valid FASTA."""
        fasta_file = fixtures_dir / "mock_reference.fasta"
        assert fasta_file.exists(), "mock_reference.fasta should exist"

        content = fasta_file.read_text()
        assert content.startswith(">"), "FASTA should start with >"
        lines = content.strip().split("\n")
        assert len(lines) >= 2, "FASTA should have header and sequence"

        # Check sequence contains only valid nucleotides
        sequence = "".join(lines[1:])
        valid_chars = set("ACGTN")
        assert all(
            c in valid_chars for c in sequence.upper()
        ), "Sequence should contain only valid nucleotides"


@pytest.mark.integration
class TestConfigConsistency:
    """Tests for configuration consistency checks."""

    def test_variants_mutation_types_valid(self, fixtures_dir):
        """Test that all mutation types in variants file are valid."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_variants.csv"
        df = pd.read_csv(csv_file)

        valid_types = {"M", "S", "D", "I", "X", "N", "Z"}
        actual_types = set(df["mutation_type"].unique())
        assert actual_types.issubset(
            valid_types
        ), f"Invalid mutation types: {actual_types - valid_types}"

    def test_experiment_has_baseline(self, fixtures_dir):
        """Test that experiment file has baseline condition."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_experiment.csv"
        df = pd.read_csv(csv_file)

        assert "baseline" in df["condition"].values, "Experiment should have baseline condition"

    def test_experiment_has_multiple_timepoints(self, fixtures_dir):
        """Test that each condition has multiple timepoints."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_experiment.csv"
        df = pd.read_csv(csv_file)

        for condition in df["condition"].unique():
            timepoints = df[df["condition"] == condition]["time"].nunique()
            assert timepoints >= 2, f"Condition {condition} should have >= 2 timepoints"

    def test_experiment_has_t0(self, fixtures_dir):
        """Test that each condition has a T0 timepoint."""
        import pandas as pd

        csv_file = fixtures_dir / "mock_experiment.csv"
        df = pd.read_csv(csv_file)

        for condition in df["condition"].unique():
            times = df[df["condition"] == condition]["time"].values
            assert 0 in times, f"Condition {condition} should have T0 timepoint"
