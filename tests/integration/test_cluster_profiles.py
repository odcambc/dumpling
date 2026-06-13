"""
Tests for the cluster execution profiles (issue #13).

Split by fidelity, per docs/cluster_profiles_design.md:
  - TestProfileConfig: the profile YAML is well-formed and selects the SLURM
    executor with the expected cross-cutting settings. Pure Python, no
    snakemake — always runs.
  - TestClusterResourcesResolve: with snakemake present, the heavy rules
    resolve their `resources: mem_mb` in a dry-run (no scheduler, no submission).
    This is the part dumpling owns; actual sbatch submission is the executor
    plugin's job and is tested upstream.
"""

import shutil
import subprocess

import pytest
import yaml


def snakemake_available():
    return shutil.which("snakemake") is not None


class TestProfileConfig:
    """The profile files exist and are valid (Tier 1 — no cluster needed)."""

    def _load(self, repo_root, name):
        path = repo_root / "workflow" / "profiles" / name / "config.yaml"
        assert path.exists(), f"Missing profile: {path}"
        with open(path) as f:
            return yaml.safe_load(f)

    def test_slurm_profile_selects_slurm_executor(self, repo_root):
        cfg = self._load(repo_root, "slurm")
        assert cfg["executor"] == "slurm"
        # Compute nodes run inside the published container.
        assert "apptainer" in cfg["software-deployment-method"]

    def test_slurm_profile_has_default_resources(self, repo_root):
        cfg = self._load(repo_root, "slurm")
        # Fallback resources for rules that don't declare their own.
        assert "mem_mb" in cfg["default-resources"]
        assert "runtime" in cfg["default-resources"]

    def test_slurm_profile_has_robustness_settings(self, repo_root):
        cfg = self._load(repo_root, "slurm")
        # Pre-emption resubmit, FS latency tolerance, mtime-only reruns.
        assert cfg["restart-times"] >= 1
        assert cfg["latency-wait"] >= 1
        assert "mtime" in cfg["rerun-triggers"]

    def test_slurm_profile_does_not_hardcode_site_values(self, repo_root):
        # account/partition are site-specific and must stay commented out, not
        # shipped with a real value (which would silently misroute jobs).
        cfg = self._load(repo_root, "slurm")
        assert "slurm_account" not in cfg["default-resources"]
        assert "slurm_partition" not in cfg["default-resources"]

    def test_default_profile_parses(self, repo_root):
        cfg = self._load(repo_root, "default")
        # Honours per-rule resources locally too (scheduled within this budget).
        assert "resources" in cfg
        assert "mtime" in cfg["rerun-triggers"]


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skipif(not snakemake_available(), reason="Snakemake not installed")
class TestClusterResourcesResolve:
    """The heavy rules resolve resources: mem_mb in a dry-run (Tier 1/2 — no
    scheduler). Validates the per-rule wiring that a SLURM run depends on."""

    def _config(self, repo_root, fixtures_dir, data_dir):
        resources_dir = repo_root / "resources"
        return yaml.safe_dump(
            {
                "experiment": "test_experiment",
                "data_dir": str(data_dir),
                "ref_dir": str(fixtures_dir),
                "experiment_file": str(fixtures_dir / "mock_experiment.csv"),
                "reference": "mock_reference.fasta",
                "variants_file": str(fixtures_dir / "mock_variants.csv"),
                "oligo_file": str(fixtures_dir / "mock_oligos.csv"),
                "orf": "1-300",
                "enrich2": False,
                "noprocess": True,
                "run_qc": False,
                "baseline_condition": "baseline",
                "remove_zeros": False,
                "regenerate_variants": False,
                "kmers": 15,
                "sam": "1.3",
                "mem": 4,
                "min_q": 30,
                "min_variant_obs": 3,
                "max_deletion_length": 3,
                "samtools_local": False,
                "rosace_local": False,
                "aligner": "bbmap",
                "adapters": str(resources_dir / "adapters.fa"),
                "contaminants": [str(resources_dir / "sequencing_artifacts.fa.gz")],
            }
        )

    def test_heavy_rules_declare_mem_mb(self, repo_root, fixtures_dir, tmp_path):
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()
        config_file = tmp_path / "config.yaml"
        config_file.write_text(self._config(repo_root, fixtures_dir, data_dir))

        # Target a sample's variant counts: pulls trim_clean_correct -> bbmap ->
        # gatk_ASM (and prepare_bbmap_index) into the DAG.
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
                "results/test_experiment/gatk/sample_A_R1_T0.variantCounts",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 0, f"dry-run failed:\n{combined}"
        # bbmap allocation (12000) and gatk allocation (6000) resolve onto the
        # jobs. Snakemake prints `mem_mb=<n>` in each job's resources line.
        assert "mem_mb=12000" in combined, f"bbmap mem_mb not resolved:\n{combined}"
        assert "mem_mb=6000" in combined, f"gatk mem_mb not resolved:\n{combined}"
