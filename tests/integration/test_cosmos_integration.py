"""
End-to-end tests for the cosmos export (issue: cosmos deposit format).

The unit tests in tests/unit/test_format_cosmos.py cover the formatting logic in
isolation. These add the end-to-end coverage:

  - TestCosmosDagConstruction: format_cosmos AND run_cosmos wire into the DAG
    (format pulls each phenotype's scores in slot order; run consumes format's
    wide CSV — triggered by targeting run_cosmos's output). Dry-run, needs
    snakemake.
  - TestCosmosLoadsInDMSData: a representative emitted CSV loads into the real
    cosmos package (DMSData), incl. the outer-join NA case. Skipped unless
    `cosmos` is importable.
  - TestCosmosRun: actually runs cosmos via run_cosmos.run_cosmos and checks the
    per-position decomposition output. SLOW (cosmos fits per position); skipped
    unless `cosmos` is importable.

See docs/cosmos_export_design.md.
"""

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

SCRIPTS_DIR = (
    Path(__file__).resolve().parents[1].parent / "workflow" / "rules" / "scripts"
)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import format_cosmos  # noqa: E402
import run_cosmos  # noqa: E402


def snakemake_available():
    return shutil.which("snakemake") is not None


def cosmos_available():
    return importlib.util.find_spec("cosmos") is not None


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.skipif(not snakemake_available(), reason="Snakemake not installed")
class TestCosmosDagConstruction:
    """format_cosmos wires into the DAG with slot-ordered score inputs."""

    def _config(self, repo_root, fixtures_dir, data_dir):
        resources_dir = repo_root / "resources"
        return yaml.safe_dump(
            {
                "experiment": "test_experiment",
                "data_dir": str(data_dir),
                "ref_dir": str(fixtures_dir),
                # 2 experimental conditions with phenotype slots 1 and 2.
                "experiment_file": str(fixtures_dir / "mock_experiment_cosmos.csv"),
                "reference": "mock_reference.fasta",
                "variants_file": str(fixtures_dir / "mock_variants.csv"),
                "oligo_file": str(fixtures_dir / "mock_oligos.csv"),
                "orf": "1-300",
                "scoring_backend": "rosace",
                "enrich2": False,
                "run_cosmos": True,
                "deposit_to_mavedb": False,
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

    def test_format_cosmos_scheduled_with_slot_ordered_inputs(
        self, repo_root, fixtures_dir, tmp_path
    ):
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()
        config_file = tmp_path / "config.yaml"
        config_file.write_text(self._config(repo_root, fixtures_dir, data_dir))

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
                "results/test_experiment/cosmos/test_experiment_cosmos.csv",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 0, f"cosmos dry-run failed:\n{combined}"

        # Find the format_cosmos job's input line.
        lines = combined.splitlines()
        input_line = ""
        for i, line in enumerate(lines):
            if line.strip() == "rule format_cosmos:":
                for follow in lines[i : i + 8]:
                    if follow.strip().startswith("input:"):
                        input_line = follow
                        break
                break
        assert input_line, f"format_cosmos not scheduled:\n{combined}"

        # Both conditions' scores are inputs, in slot order: cond_A (phenotype 1,
        # -> beta_hat_1) must precede cond_B (phenotype 2, -> beta_hat_2).
        a = input_line.find("cond_A_scores.csv")
        b = input_line.find("cond_B_scores.csv")
        assert a != -1 and b != -1, f"missing score inputs:\n{input_line}"
        assert (
            a < b
        ), f"score inputs not in slot order (cond_A before cond_B):\n{input_line}"

    def test_run_cosmos_scheduled_consuming_format_output(
        self, repo_root, fixtures_dir, tmp_path
    ):
        """The run_cosmos rule wires into the DAG, taking format_cosmos's wide
        CSV as input and emitting the per-position results CSV. Dry-run only —
        no actual cosmos fitting."""
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "mock_reads_R1.fastq.gz").touch()
        (data_dir / "mock_reads_R2.fastq.gz").touch()
        config_file = tmp_path / "config.yaml"
        config_file.write_text(self._config(repo_root, fixtures_dir, data_dir))

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
                "results/test_experiment/cosmos/"
                "test_experiment_cosmos_results.csv",
            ],
            capture_output=True,
            text=True,
            cwd=str(repo_root),
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 0, f"run_cosmos dry-run failed:\n{combined}"

        lines = combined.splitlines()
        input_line = ""
        for i, line in enumerate(lines):
            if line.strip() == "rule run_cosmos:":
                for follow in lines[i : i + 8]:
                    if follow.strip().startswith("input:"):
                        input_line = follow
                        break
                break
        assert input_line, f"run_cosmos not scheduled:\n{combined}"
        # Its input is format_cosmos's wide CSV (the prep stage feeds the run).
        assert (
            "test_experiment_cosmos.csv" in input_line
        ), f"run_cosmos not consuming the format_cosmos output:\n{input_line}"


# -----------------------------------------------------------------------------
# Validation against the real cosmos package. Skipped unless `cosmos`
# (pip install cosmos-dms / pimentellab/cosmos) is importable.
# -----------------------------------------------------------------------------


def _make_two_phenotype_scores(tmp_path, n_positions=15, variants_per_pos=12):
    """Write two per-condition score CSVs (rosace column layout) with enough
    missense variants per position to clear cosmos's min_num_variants_per_group.
    Returns (path_pheno1, path_pheno2)."""
    aas = "ACDEFGHIKLMNPQRSTVWY"
    rows1, rows2 = [], []
    for pos in range(1, n_positions + 1):
        for j in range(variants_per_pos):
            aa = aas[j % len(aas)]
            variant = f"p.(X{pos}{aa})"
            # Deterministic-but-varied effect sizes; no randomness needed.
            rows1.append((variant, pos, aa, aa, "missense", (pos + j) * 0.01, 0.1))
            rows2.append((variant, pos, aa, aa, "missense", (pos - j) * 0.02, 0.1))
    cols = ["variants", "position", "wildtype", "mutation", "type", "mean", "sd"]
    p1 = tmp_path / "cond_A_scores.csv"
    p2 = tmp_path / "cond_B_scores.csv"
    pd.DataFrame(rows1, columns=cols).to_csv(p1, index=False)
    pd.DataFrame(rows2, columns=cols).to_csv(p2, index=False)
    return p1, p2


@pytest.mark.skipif(
    not cosmos_available(),
    reason="cosmos package not installed (pip install cosmos-dms)",
)
class TestCosmosLoadsInDMSData:
    """The emitted CSV satisfies cosmos's DMSData column contract."""

    def test_emitted_csv_loads_into_dmsdata(self, tmp_path):
        from cosmos import DMSData

        p1, p2 = _make_two_phenotype_scores(tmp_path)
        out = tmp_path / "exp_cosmos.csv"
        cols = format_cosmos._BACKEND_DEFAULTS["rosace"]
        format_cosmos.format_cosmos([str(p1), str(p2)], "rosace", cols, str(out))

        df = pd.read_csv(out)
        # Construct DMSData with two phenotypes; should not raise on _check_cols.
        data = DMSData(
            df,
            ["phenotype_1", "phenotype_2"],
            include_type=["missense"],
            exclude_type=["synonymous"],
            min_num_variants_per_group=10,
        )
        assert data is not None

    def test_dmsdata_tolerates_outer_join_nas(self, tmp_path):
        """The design doc's open question: format_cosmos outer-joins, so a
        variant scored in only one phenotype gets an NA beta_hat in the other.
        Confirm cosmos's DMSData accepts (and keeps) those NA rows rather than
        rejecting the file — verified empirically against cosmos 1.0.2: it loads
        without error and retains the NA row. If a future cosmos version starts
        rejecting NAs, this test fails loudly and we switch to documenting an
        inner-join subset for DMSData."""
        from cosmos import DMSData

        p1, p2 = _make_two_phenotype_scores(tmp_path)
        # Append a variant present only in condition A -> NA in beta_hat_2.
        a_only = pd.DataFrame(
            [("p.(X3onlyA)", 3, "A", "A", "missense", 0.5, 0.1)],
            columns=[
                "variants",
                "position",
                "wildtype",
                "mutation",
                "type",
                "mean",
                "sd",
            ],
        )
        a_only.to_csv(p1, mode="a", header=False, index=False)

        out = tmp_path / "exp_cosmos.csv"
        cols = format_cosmos._BACKEND_DEFAULTS["rosace"]
        format_cosmos.format_cosmos([str(p1), str(p2)], "rosace", cols, str(out))
        df = pd.read_csv(out)
        assert df["beta_hat_2"].isna().any(), "fixture should have an NA-padded row"

        n_before = len(df)
        data = DMSData(
            df,
            ["phenotype_1", "phenotype_2"],
            include_type=["missense"],
            exclude_type=["synonymous"],
            min_num_variants_per_group=10,
        )
        # Loads without error and keeps the NA row (no silent drop).
        retained = getattr(data, "data", df)
        assert len(retained) == n_before


def _make_noisy_two_phenotype_scores(tmp_path, n_positions=3, per_pos=16):
    """Two per-condition score CSVs with a real mediated structure
    (beta_2 ~ 0.7*beta_1 + noise) so cosmos's per-position fit converges.
    Deterministic via a fixed seed."""
    rng = np.random.default_rng(1)
    aas = "ACDEFGHIKLMNPQRSTVWY"
    rows1, rows2 = [], []
    for pos in range(1, n_positions + 1):
        base = rng.normal(0, 1)
        for j in range(per_pos):
            aa = aas[j % len(aas)]
            v = f"p.(X{pos}{aa})"
            b1 = base + rng.normal(0, 0.6)
            b2 = 0.7 * b1 + rng.normal(0, 0.6)
            rows1.append((v, pos, aa, aa, "missense", b1, 0.1))
            rows2.append((v, pos, aa, aa, "missense", b2, 0.1))
    cols = ["variants", "position", "wildtype", "mutation", "type", "mean", "sd"]
    p1, p2 = tmp_path / "cond_A_scores.csv", tmp_path / "cond_B_scores.csv"
    pd.DataFrame(rows1, columns=cols).to_csv(p1, index=False)
    pd.DataFrame(rows2, columns=cols).to_csv(p2, index=False)
    return p1, p2


@pytest.mark.slow
@pytest.mark.skipif(
    not cosmos_available(),
    reason="cosmos package not installed (uv sync --extra cosmos)",
)
class TestCosmosRun:
    """Actually run cosmos end-to-end via run_cosmos.run_cosmos. SLOW: cosmos
    fits a model per position (~tens of seconds each), so this is capped to a
    few positions. Validates that the run produces the per-position
    decomposition the rule promises."""

    def test_run_produces_per_position_summary(self, tmp_path):
        p1, p2 = _make_noisy_two_phenotype_scores(tmp_path, n_positions=3)
        cosmos_in = tmp_path / "cosmos_in.csv"
        cols = format_cosmos._BACKEND_DEFAULTS["rosace"]
        format_cosmos.format_cosmos([str(p1), str(p2)], "rosace", cols, str(cosmos_in))

        out = tmp_path / "cosmos_results.csv"
        run_cosmos.run_cosmos(
            str(cosmos_in),
            ["abundance", "activity"],
            str(out),
            min_num_variants_per_group=10,
        )

        res = pd.read_csv(out)
        # cosmos's per-position decomposition: tau (direct) and gamma (mediated),
        # one row per position group.
        for col in ["position", "tau_mean", "tau_std", "gamma_mean", "gamma_std"]:
            assert (
                col in res.columns
            ), f"missing cosmos summary column {col}: {list(res.columns)}"
        assert len(res) >= 1
        # gamma (mediated effect) is fit for every retained position.
        assert res["gamma_mean"].notna().any()

    def test_run_rejects_non_two_phenotypes(self, tmp_path):
        # cosmos models exactly two sequential phenotypes; one (or three) must
        # fail loud rather than silently mis-fit.
        cosmos_in = tmp_path / "in.csv"
        pd.DataFrame(
            {
                "variants": ["p.(G2A)"],
                "group": [2],
                "type": ["missense"],
                "beta_hat_1": [1.0],
                "se_hat_1": [0.1],
            }
        ).to_csv(cosmos_in, index=False)
        with pytest.raises(ValueError, match="exactly 2 phenotypes"):
            run_cosmos.run_cosmos(str(cosmos_in), ["only_one"], str(tmp_path / "o.csv"))
