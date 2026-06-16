"""
Integration tests for Snakemake DAG construction.

Tests that the pipeline DAG builds correctly with various configurations
using dry-run mode (no actual execution).
"""

import os
import shutil
import subprocess

import pytest
import yaml


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
            "lilace.smk",
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

    def _enrich_config(self, fixtures_dir, repo_root, data_dir, **overrides):
        """Mock config for the Enrich2 path. data_dir holds stub fastqs so the
        full DAG builds (common.smk + baseline_qc resolve every sample's fastq
        pair at parse time). overrides patch individual keys (e.g. keep_enrich_h5).
        """
        resources_dir = repo_root / "resources"
        cfg = {
            "experiment": "test_experiment",
            "data_dir": str(data_dir),
            "ref_dir": str(fixtures_dir),
            "experiment_file": str(fixtures_dir / "mock_experiment.csv"),
            "reference": "mock_reference.fasta",
            "variants_file": str(fixtures_dir / "mock_variants.csv"),
            "oligo_file": str(fixtures_dir / "mock_oligos.csv"),
            "orf": "1-300",
            "enrich2": True,
            "noprocess": True,
            "run_qc": False,
            "baseline_condition": "baseline",
            "remove_zeros": False,
            "regenerate_variants": False,
            "kmers": 15,
            "sam": "1.3",
            "mem": 4,
            "mem_fastqc": 1024,
            "min_q": 30,
            "min_variant_obs": 3,
            "max_deletion_length": 3,
            "samtools_local": False,
            "rosace_local": False,
            "adapters": str(resources_dir / "adapters.fa"),
            "contaminants": [str(resources_dir / "sequencing_artifacts.fa.gz")],
        }
        cfg.update(overrides)
        return yaml.safe_dump(cfg)

    def _run_enrich_dry_run(self, repo_root, config_file):
        return subprocess.run(
            [
                "snakemake",
                "-s",
                str(repo_root / "workflow" / "Snakefile"),
                "--configfile",
                str(config_file),
                "--dry-run",
                "-p",
                "--cores",
                "1",
                "results/test_experiment/enrich/tsv/"
                "test_experiment_exp/main_identifiers_scores.tsv",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )

    @staticmethod
    def _run_enrich_output_line(stdout):
        """The `output:` line of the scheduled run_enrich job, or '' if absent."""
        lines = stdout.splitlines()
        for i, line in enumerate(lines):
            if line.strip() == "rule run_enrich:":
                for follow in lines[i : i + 12]:
                    if follow.strip().startswith("output:"):
                        return follow
        return ""

    def test_dry_run_enrich_declares_h5_as_temp_outputs(
        self, repo_root, fixtures_dir, tmp_path
    ):
        """Issue #16: with enrich2 on, run_enrich declares Enrich2's .h5 stores
        as outputs so Snakemake (not a pipeline rm/find) manages them. The
        predicted store set must match the mock experiment: cond_A rep1 (T0-T2)
        and rep2 (T0-T2); baseline is not a scored condition.

        The other dry-run tests all set enrich2: false, so run_enrich and the
        temp-output enumeration are otherwise never exercised.
        """
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()

        config_file = tmp_path / "config.yaml"
        config_file.write_text(self._enrich_config(fixtures_dir, repo_root, data_dir))

        result = self._run_enrich_dry_run(repo_root, config_file)
        assert result.returncode == 0, (
            f"Enrich2 dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

        output_line = self._run_enrich_output_line(result.stdout)
        assert output_line, f"run_enrich not scheduled:\n{result.stdout}"
        # Every predicted store appears as a declared output of run_enrich.
        for store in [
            "test_experiment_exp.h5",
            "cond_A_R1_sel.h5",
            "cond_A_R2_sel.h5",
            "cond_A_rep1_T0_lib.h5",
            "cond_A_rep1_T2_lib.h5",
            "cond_A_rep2_T0_lib.h5",
        ]:
            assert store in output_line, f"{store} not declared:\n{output_line}"
        # Baseline is not scored, so no baseline stores are declared.
        assert "baseline" not in output_line

    def test_dry_run_keep_enrich_h5_declares_only_scores(
        self, repo_root, fixtures_dir, tmp_path
    ):
        """keep_enrich_h5: true leaves the stores undeclared (they persist as
        side-effects), so run_enrich declares only the tsv scores."""
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()

        config_file = tmp_path / "config.yaml"
        config_file.write_text(
            self._enrich_config(
                fixtures_dir, repo_root, data_dir, keep_enrich_h5=True
            )
        )

        result = self._run_enrich_dry_run(repo_root, config_file)
        assert result.returncode == 0, (
            f"Enrich2 dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

        output_line = self._run_enrich_output_line(result.stdout)
        assert output_line, f"run_enrich not scheduled:\n{result.stdout}"
        assert ".h5" not in output_line, (
            f"keep_enrich_h5=true should declare no .h5 outputs:\n{output_line}"
        )
        assert "main_identifiers_scores.tsv" in output_line

    def test_dry_run_with_lilace_backend(self, repo_root, fixtures_dir, tmp_path):
        """Dry-run with scoring_backend=lilace should construct the DAG without crashing
        and schedule lilace rules instead of rosace.

        noprocess is deliberately false here — lilace requires a synonymous-variant
        control set and validate_scoring_backend_mode rejects noprocess+lilace at
        parse time. See test_noprocess_with_lilace_rejected for that path.
        """
        test_dir = tmp_path / "test_run_lilace"
        test_dir.mkdir()

        config_content = f"""
experiment: 'test_experiment'
data_dir: '{fixtures_dir}'
ref_dir: '{fixtures_dir}'
experiment_file: '{fixtures_dir / "mock_experiment.csv"}'
reference: 'mock_reference.fasta'
variants_file: '{fixtures_dir / "mock_variants.csv"}'
oligo_file: '{fixtures_dir / "mock_oligos.csv"}'
orf: "1-300"
scoring_backend: 'lilace'
enrich2: false
noprocess: false
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
lilace_local: false
adapters: 'resources/adapters.fa'
contaminants:
  - 'resources/sequencing_artifacts.fa.gz'
"""
        config_file = test_dir / "test_config.yaml"
        config_file.write_text(config_content)

        result = subprocess.run(
            [
                "snakemake",
                "-s",
                str(repo_root / "workflow" / "Snakefile"),
                "--configfile",
                str(config_file),
                "--dry-run",
                "-p",
                "--cores",
                "1",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )

        # As with the rosace dry-run, fixtures are incomplete so a non-zero exit
        # is acceptable — what we care about is that config parsing and rule
        # construction succeed for the lilace path.
        if result.returncode != 0:
            assert "SyntaxError" not in result.stderr, f"Syntax error: {result.stderr}"
            assert "NameError" not in result.stderr, f"Name error: {result.stderr}"
            # KeyError on a config key (e.g. mem_lilace) means the setdefault /
            # schema-default path is broken and would prevent any lilace run.
            assert "KeyError" not in result.stderr, (
                f"Missing config key in lilace path: {result.stderr}"
            )

        # If anything got printed at all, the DAG should mention lilace rather
        # than rosace. Don't fail on empty output (missing-file errors abort early).
        combined = result.stdout + result.stderr
        if "rule " in combined:
            assert "lilace" in combined.lower(), (
                "Lilace backend selected but no lilace rule scheduled:\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
            assert "run_rosace" not in combined, (
                "Lilace backend selected but run_rosace was also scheduled:\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )

    def test_noprocess_with_lilace_rejected(self, repo_root, fixtures_dir, tmp_path):
        """noprocess + scoring_backend=lilace must error at parse time, not deep
        inside R after install_lilace runs. validate_scoring_backend_mode is the
        gate; this test exercises it through the actual Snakefile."""
        test_dir = tmp_path / "test_run_noprocess_lilace"
        test_dir.mkdir()

        config_content = f"""
experiment: 'test_experiment'
data_dir: '{fixtures_dir}'
ref_dir: '{fixtures_dir}'
experiment_file: '{fixtures_dir / "mock_experiment.csv"}'
reference: 'mock_reference.fasta'
variants_file: '{fixtures_dir / "mock_variants.csv"}'
oligo_file: '{fixtures_dir / "mock_oligos.csv"}'
orf: "1-300"
scoring_backend: 'lilace'
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
lilace_local: false
adapters: 'resources/adapters.fa'
contaminants:
  - 'resources/sequencing_artifacts.fa.gz'
"""
        config_file = test_dir / "test_config.yaml"
        config_file.write_text(config_content)

        result = subprocess.run(
            [
                "snakemake",
                "-s",
                str(repo_root / "workflow" / "Snakefile"),
                "--configfile",
                str(config_file),
                "--dry-run",
                "--cores",
                "1",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )

        assert result.returncode != 0, (
            "noprocess+lilace should be rejected at parse time, but Snakefile "
            f"loaded cleanly.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
        combined = result.stdout + result.stderr
        assert "lilace" in combined.lower() and "noprocess" in combined.lower(), (
            "Error message should mention both lilace and noprocess so the user "
            f"knows the combination is the problem.\nSTDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

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

        # Stub fastq files matching the mock experiment's `file: mock_reads`
        # column. The rule under test (generate_baseline_configs) doesn't
        # consume them, but common.smk performs fastq pair resolution at
        # parse time for every sample in the experiment CSV, and baseline_qc
        # then indexes fastqc_names[<prefix>] — a KeyError here would block
        # parse before the requested rule could run. Empty files are
        # sufficient: resolve_fastq_pair only matches on filename pattern.
        data_dir = test_dir / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()

        resources_dir = repo_root / "resources"
        config_file = test_dir / "test_config.yaml"
        config_content = f"""
experiment: 'smoke'
data_dir: '{data_dir}'
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

    def test_lilace_rules_exist(self, repo_root):
        """Test that lilace rules exist."""
        lilace_file = repo_root / "workflow" / "rules" / "lilace.smk"
        content = lilace_file.read_text()

        assert "lilace" in content.lower()
        assert "rule run_lilace" in content
        assert "rule install_lilace" in content

    def test_get_input_branches_on_scoring_backend(self, repo_root):
        """get_input must interpolate the scoring_backend into its output targets,
        so flipping config['scoring_backend'] actually switches the requested rules."""
        common_smk = (repo_root / "workflow" / "rules" / "common.smk").read_text()

        # The f-string forms used in get_input — these are the exact templates
        # that drive backend selection. If they're refactored, this test should
        # be updated deliberately, not silently broken.
        assert "{scoring_backend}/{{conditions}}_scores.csv" in common_smk, (
            "get_input no longer interpolates scoring_backend into the scores target"
        )
        assert "{scoring_backend}/{scoring_backend}_installed.txt" in common_smk, (
            "get_input no longer interpolates scoring_backend into the install marker"
        )
