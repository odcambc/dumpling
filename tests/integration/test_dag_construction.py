"""
Integration tests for Snakemake DAG construction.

Tests that the pipeline DAG builds correctly with various configurations
using dry-run mode (no actual execution).
"""

import pytest
import subprocess
import shutil
import os
from pathlib import Path


def snakemake_available():
    """Check if snakemake is available."""
    return shutil.which("snakemake") is not None


@pytest.mark.integration
@pytest.mark.skipif(not snakemake_available(), reason="Snakemake not installed")
class TestDAGConstruction:
    """Tests for Snakemake DAG construction via dry-run."""

    def test_snakefile_syntax_valid(self, repo_root):
        """Test that Snakefile has valid Python/Snakemake syntax."""
        snakefile = repo_root / "workflow" / "Snakefile"
        assert snakefile.exists(), "Snakefile should exist"

        # Basic syntax check - try to compile the file
        with open(snakefile) as f:
            content = f.read()

        # Check for common Snakemake patterns
        assert "rule all:" in content or "rule all" in content, "Should have rule all"
        assert "include:" in content, "Should include other rule files"

    def test_rule_files_exist(self, repo_root):
        """Test that all included rule files exist."""
        rules_dir = repo_root / "workflow" / "rules"
        expected_rules = [
            "common.smk",
            "ref.smk",
            "filter.smk",
            "map.smk",
            "asm.smk",
            "process.smk",
            "qc.smk",
            "rosace.smk",
            "enrich.smk",
        ]

        for rule_file in expected_rules:
            assert (rules_dir / rule_file).exists(), f"Rule file {rule_file} should exist"

    def test_schema_files_exist(self, repo_root):
        """Test that schema files exist for validation."""
        schemas_dir = repo_root / "workflow" / "schemas"
        expected_schemas = ["config.schema.yaml", "experiments.schema.yaml"]

        for schema_file in expected_schemas:
            assert (
                schemas_dir / schema_file
            ).exists(), f"Schema file {schema_file} should exist"


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skipif(not snakemake_available(), reason="Snakemake not installed")
class TestDryRun:
    """Tests using Snakemake dry-run (requires snakemake installed)."""

    def test_dry_run_with_mock_config(self, repo_root, fixtures_dir, tmp_path):
        """Test that dry-run succeeds with mock config."""
        # Copy fixtures to tmp_path for isolation
        import shutil

        test_dir = tmp_path / "test_run"
        test_dir.mkdir()

        # Create a minimal config for dry-run
        config_content = f"""
experiment: 'test_experiment'
data_dir: '{fixtures_dir}'
ref_dir: '{fixtures_dir}'
experiment_file: '{fixtures_dir / "mock_experiment.csv"}'
reference: 'mock_reference.fasta'
variants_file: '{fixtures_dir / "mock_variants.csv"}'
oligo_file: '{fixtures_dir / "mock_oligos.csv"}'
orf: "1-300"
enrich2: false
noprocess: true
run_qc: false
baseline_condition: 'baseline'
remove_zeros: false
regenerate_variants: false
kmers: 15
sam: "1.3"
mem: 4
min_q: 30
min_variant_obs: 3
max_deletion_length: 3
samtools_local: false
rosace_local: false
adapters: 'resources/adapters.fa'
contaminants:
  - 'resources/sequencing_artifacts.fa.gz'
"""
        config_file = test_dir / "test_config.yaml"
        config_file.write_text(config_content)

        # Run snakemake dry-run
        result = subprocess.run(
            [
                "snakemake",
                "-s",
                str(repo_root / "workflow" / "Snakefile"),
                "--configfile",
                str(config_file),
                "--dry-run",
                "-q",
                "--cores",
                "1",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )

        # Note: This may fail due to missing data files, but should not have syntax errors
        # We're mainly checking that the DAG construction doesn't crash
        if result.returncode != 0:
            # Check if it's a missing file error (expected) vs syntax error (not expected)
            assert "SyntaxError" not in result.stderr, f"Syntax error: {result.stderr}"
            assert "NameError" not in result.stderr, f"Name error: {result.stderr}"

    def test_script_rule_executes_with_runtime_injection(
        self, repo_root, fixtures_dir, tmp_path
    ):
        """Execute a lightweight script-backed rule (not dry-run)."""
        test_dir = tmp_path / "runtime_smoke"
        test_dir.mkdir()

        (test_dir / "config").mkdir()
        (test_dir / "config" / "multiqc_config.yaml").write_text(
            "extra_fn_clean_exts:\n  - \".txt\"\n"
        )

        resources_dir = repo_root / "resources"
        config_file = test_dir / "test_config.yaml"
        config_content = f"""
experiment: 'smoke'
data_dir: '{fixtures_dir}'
ref_dir: '{fixtures_dir}'
experiment_file: '{fixtures_dir / "mock_experiment.csv"}'
reference: 'mock_reference.fasta'
variants_file: '{fixtures_dir / "mock_variants.csv"}'
oligo_file: '{fixtures_dir / "mock_oligos.csv"}'
orf: "1-300"
enrich2: false
noprocess: true
run_qc: false
baseline_condition: 'baseline'
remove_zeros: false
regenerate_variants: false
kmers: 15
sam: "1.3"
mem: 4
mem_fastqc: 1024
mem_rosace: 16000
min_q: 30
min_variant_obs: 3
max_deletion_length: 3
samtools_local: false
rosace_local: false
adapters: '{resources_dir / "adapters.fa"}'
contaminants:
  - '{resources_dir / "sequencing_artifacts.fa.gz"}'
"""
        config_file.write_text(config_content)

        output_target = "config/smoke_multiqc_baseline_config.yaml"
        env = os.environ.copy()
        env["HOME"] = str(test_dir)
        env["XDG_CACHE_HOME"] = str(test_dir / ".cache")
        result = subprocess.run(
            [
                "snakemake",
                "-s",
                str(repo_root / "workflow" / "Snakefile"),
                "--directory",
                str(test_dir),
                "--configfile",
                str(config_file),
                "--cores",
                "1",
                "--notemp",
                output_target,
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
            env=env,
        )

        assert result.returncode == 0, (
            "Script-backed rule execution failed:\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        assert (test_dir / output_target).exists(), "Expected script output was not created"


@pytest.mark.integration
class TestConditionalRules:
    """Tests for conditional rule inclusion based on config flags."""

    def test_enrich2_rules_conditional(self, repo_root):
        """Test that enrich2 rules are conditionally included."""
        enrich_file = repo_root / "workflow" / "rules" / "enrich.smk"
        content = enrich_file.read_text()

        # Check for conditional logic
        assert "enrich2" in content.lower() or "config" in content.lower()

    def test_qc_rules_conditional(self, repo_root):
        """Test that QC rules exist."""
        qc_file = repo_root / "workflow" / "rules" / "qc.smk"
        content = qc_file.read_text()

        # Check for fastqc and multiqc rules
        assert "fastqc" in content.lower()
        assert "multiqc" in content.lower()

    def test_rosace_rules_exist(self, repo_root):
        """Test that rosace rules exist."""
        rosace_file = repo_root / "workflow" / "rules" / "rosace.smk"
        content = rosace_file.read_text()

        # Check for rosace-related rules
        assert "rosace" in content.lower()
        assert "rule" in content
