"""
Shared pytest fixtures for Dumpling pipeline tests.

This module provides common fixtures for unit and integration tests including:
- Mock Snakemake objects
- Sample DataFrames (designed variants, experiments)
- Reference sequences and GATK outputs
- Temporary directory management
"""

import os
import sys
import pytest
import numpy
import pandas as pd
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock

# Add workflow scripts to path for imports
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "workflow" / "rules" / "scripts"))


def _make_bootstrap_snakemake():
    """Create a minimal snakemake-like object for module import time."""
    fixtures_dir = REPO_ROOT / "tests" / "fixtures"
    bootstrap_ref = REPO_ROOT / "references" / "cftr_d508_frag1.fasta"
    bootstrap_oligos = REPO_ROOT / "config" / "oligos" / "CFTR_1_d508_oligos.csv"
    return MagicMock(
        config={
            "experiment": "mock_experiment",
            "experiment_file": str(fixtures_dir / "mock_experiment.csv"),
            "baseline_condition": "baseline",
            "variants_file": str(fixtures_dir / "mock_variants.csv"),
            "orf": "229-2817",
        },
        input={
            "oligo_file": str(bootstrap_oligos),
            "ref_fasta": str(bootstrap_ref),
        },
        output=[os.devnull],
        params={},
        log=[os.devnull],
    )


def _install_snakemake_import_shim():
    """Ensure ``from snakemake.script import snakemake`` works during pytest imports."""
    bootstrap = _make_bootstrap_snakemake()

    try:
        import snakemake.script as snakemake_script

        if not hasattr(snakemake_script, "snakemake"):
            snakemake_script.snakemake = bootstrap
        return
    except ModuleNotFoundError:
        pass

    snakemake_module = ModuleType("snakemake")
    script_module = ModuleType("snakemake.script")
    script_module.snakemake = bootstrap
    snakemake_module.script = script_module
    sys.modules.setdefault("snakemake", snakemake_module)
    sys.modules["snakemake.script"] = script_module


_install_snakemake_import_shim()


# =============================================================================
# Mock Snakemake Object
# =============================================================================


class MockSnakemake:
    """
    Mock Snakemake object for testing scripts that depend on snakemake context.

    Usage:
        snakemake = MockSnakemake(
            config={'experiment': 'test'},
            input={'ref': 'path/to/ref.fasta'},
            output=['output.csv'],
            log=['/dev/null']
        )
    """

    def __init__(
        self,
        config=None,
        input=None,
        output=None,
        params=None,
        log=None,
        wildcards=None,
        threads=1,
        resources=None,
    ):
        self.config = config or {}
        self.input = self._to_namespace(input or {})
        self.output = output or []
        self.params = self._to_namespace(params or {})
        self.log = log or ["/dev/null"]
        self.wildcards = self._to_namespace(wildcards or {})
        self.threads = threads
        self.resources = self._to_namespace(resources or {})

    def _to_namespace(self, d):
        """Convert dict to object with attribute access."""
        if isinstance(d, dict):
            ns = MagicMock()
            for k, v in d.items():
                setattr(ns, k, v)
            # Allow dict-style access too
            ns.__getitem__ = lambda self, key: getattr(self, key)
            ns.get = lambda key, default=None: getattr(ns, key, default)
            return ns
        return d


@pytest.fixture
def mock_snakemake():
    """Factory fixture for creating mock Snakemake objects."""
    return MockSnakemake


# =============================================================================
# Path Fixtures
# =============================================================================


@pytest.fixture
def fixtures_dir():
    """Path to the test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def repo_root():
    """Path to the repository root."""
    return REPO_ROOT


# =============================================================================
# Reference Sequence Fixtures
# =============================================================================


@pytest.fixture
def mock_reference_sequence():
    """
    A synthetic 500bp DNA reference sequence for testing.
    Contains a clear ORF from position 1-300 (100 codons).
    """
    # Start with ATG, end with TAA, varied codons in between
    sequence = (
        "ATGGCTAAGGTCGCTTTCTGGCTGCTGCTGTCGAAGCGTGACAACCTGACC"  # 50bp
        "TACGAGTACAAGAACGGCTACCTGAACGCTCGTGACGTCAACAAGCTGTAC"  # 100bp
        "GACTTCAAGAACGCTGACCTGAACGGTCGTAACAAGCTGAACGACTTCAAG"  # 150bp
        "AACGCTGACCTGAACGGTCGTAACAAGCTGAACGACTTCAAGAACGCTGAC"  # 200bp
        "CTGAACGGTCGTAACAAGCTGAACGACTTCAAGAACGCTGACCTGAACGGT"  # 250bp
        "CGTAACAAGCTGAACGACTTCAAGAACGCTGACCTGAACGGTCGTAACTAA"  # 300bp (stop)
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"  # 350bp (UTR)
        "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"  # 400bp
        "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"  # 450bp
        "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"  # 500bp
    )
    return sequence


@pytest.fixture
def mock_aa_sequence():
    """
    Amino acid sequence corresponding to mock_reference_sequence ORF (codons 1-99).
    """
    return "MAKVAFWLLLLSKRDNLTYEYKNGYLAARDVNKLYDKFNADLNGRNNKLNDFKNADLNGRNNKLNDFKNADLNGRNNKLNDFKNADLNGRNNKLNDFKNADLNGR"


@pytest.fixture
def short_ref_aa_sequence():
    """Short reference AA sequence for simple tests."""
    return "MKVAFWLLLS"


# =============================================================================
# Designed Variants Fixtures
# =============================================================================


@pytest.fixture
def designed_variants_df():
    """
    Small mock DataFrame representing designed variants.
    Each row corresponds to a single variant the pipeline expects.
    """
    data = {
        "count": [0, 0, 0, 0, 0],
        "pos": [10, 11, 12, 13, 14],
        "mutation_type": ["M", "S", "D", "I", "X"],
        "name": ["G10A", "R11R", "P12del", "L13_K14insA", "W14X"],
        "codon": ["GCT", "CGT", "", "GCT", "TGA"],
        "mutation": ["A", "R", "D_1", "I_1", "X"],
        "length": [1, 1, 1, 1, 1],
        "hgvs": ["p.(G10A)", "p.(R11R)", "p.(P12del)", "p.(L13_K14insA)", "p.(W14X)"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def full_designed_variants_df():
    """
    Larger mock DataFrame with multiple mutation types at each position.
    More realistic for integration testing.
    """
    variants = []
    for pos in range(1, 21):
        # Missense variants
        for aa in ["A", "C", "D", "E", "F"]:
            variants.append({
                "count": 0,
                "pos": pos,
                "mutation_type": "M",
                "name": f"X{pos}{aa}",
                "codon": "GCT",
                "mutation": aa,
                "length": 1,
                "hgvs": f"p.(X{pos}{aa})",
            })
        # Synonymous
        variants.append({
            "count": 0,
            "pos": pos,
            "mutation_type": "S",
            "name": f"X{pos}X",
            "codon": "XXX",
            "mutation": "S",
            "length": 1,
            "hgvs": f"p.(X{pos}=)",
        })
        # Deletion
        variants.append({
            "count": 0,
            "pos": pos,
            "mutation_type": "D",
            "name": f"X{pos}del",
            "codon": "",
            "mutation": "D_1",
            "length": 1,
            "hgvs": f"p.(X{pos}del)",
        })
    return pd.DataFrame(variants)


# =============================================================================
# GATK Output Fixtures
# =============================================================================


@pytest.fixture
def mock_gatk_list():
    """
    Mock GATK ASM output as a list of lists (parsed variantCounts format).
    Columns: counts, coverage, mean_length, length_NT, NT, length_codon, codon, AA, mutations
    """
    return [
        ["5", "100", "0.1", "1", "1774:A>C", "1", "10:CAA>GCT", "M:G>A", "G10A"],
        ["2", "100", "1.2", "2", "1774:A>C,1823:G>A", "1", "11:CGC>CGT", "S:R>R", ""],
        ["10", "100", "0", "3", "1774:A>-,1775:A>-,1776:A>-", "1", "12:AAA>---", "D:P>-", "P12del"],
        ["3", "100", "0.5", "1", "1780:G>A", "1", "14:TGG>TGA", "N:W>*", "W14X"],
    ]


@pytest.fixture
def mock_gatk_file(tmp_path):
    """Create a temporary GATK variantCounts file."""
    gatk_content = """5	100	0.1	1	1774:A>C	1	10:CAA>GCT	M:G>A	G10A
2	100	1.2	2	1774:A>C,1823:G>A	1	11:CGC>CGT	S:R>R
10	100	0	3	1774:A>-,1775:A>-,1776:A>-	1	12:AAA>---	D:P>-	P12del
3	100	0.5	1	1780:G>A	1	14:TGG>TGA	N:W>*	W14X
7	100	0.3	1	1785:C>T	1	15:AAA>AAG	M:K>K
"""
    gatk_file = tmp_path / "test.variantCounts"
    gatk_file.write_text(gatk_content)
    return gatk_file


# =============================================================================
# Experiment Configuration Fixtures
# =============================================================================


@pytest.fixture
def mock_experiment_df():
    """
    Mock experiment DataFrame matching the experiments.schema.yaml format.
    """
    data = {
        "sample": ["sample_A_T0", "sample_A_T1", "sample_B_T0", "sample_B_T1"],
        "condition": ["cond_A", "cond_A", "baseline", "baseline"],
        "replicate": [1, 1, 1, 1],
        "time": [0, 1, 0, 1],
        "tile": [1, 1, 1, 1],
        "file": ["mock_reads", "mock_reads", "mock_reads", "mock_reads"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def mock_config():
    """
    Mock configuration dictionary matching config.schema.yaml.
    """
    return {
        "experiment": "test_experiment",
        "experiment_file": "config/mock_experiment.csv",
        "data_dir": "tests/fixtures",
        "ref_dir": "tests/fixtures",
        "reference": "mock_reference.fasta",
        "variants_file": "tests/fixtures/mock_variants.csv",
        "oligo_file": "tests/fixtures/mock_oligos.csv",
        "orf": "1-300",
        "enrich2": True,
        "remove_zeros": True,
        "regenerate_variants": False,
        "noprocess": False,
        "run_qc": False,
        "baseline_condition": "baseline",
        "max_deletion_length": 3,
        "kmers": 15,
        "sam": "1.3",
        "mem": 16,
        "min_q": 30,
        "min_variant_obs": 3,
        "samtools_local": False,
        "rosace_local": False,
        "bbtools_use_bgzip": True,
        "adapters": "resources/adapters.fa",
        "contaminants": [
            "resources/sequencing_artifacts.fa.gz",
            "resources/phix174_ill.ref.fa.gz",
        ],
    }


# =============================================================================
# Enrich Format Fixtures
# =============================================================================


@pytest.fixture
def mock_enrich_files(tmp_path):
    """
    Create mock Enrich2 format TSV files for testing remove_zeros.
    """
    data_dir = tmp_path / "enrich_format"
    data_dir.mkdir()

    # Sample 1: varA=10, varB=0, varC=5
    sample1 = data_dir / "sample1.tsv"
    sample1.write_text("hgvs\tcount\nvarA\t10\nvarB\t0\nvarC\t5\n")

    # Sample 2: varA=0, varB=0, varC=5, varD=0
    sample2 = data_dir / "sample2.tsv"
    sample2.write_text("hgvs\tcount\nvarA\t0\nvarB\t0\nvarC\t5\nvarD\t0\n")

    return data_dir


# =============================================================================
# File System Fixtures
# =============================================================================


@pytest.fixture
def mock_reference_fasta(tmp_path, mock_reference_sequence):
    """Create a temporary FASTA reference file."""
    fasta_file = tmp_path / "mock_reference.fasta"
    fasta_file.write_text(f">mock_reference\n{mock_reference_sequence}\n")
    return fasta_file


@pytest.fixture
def mock_variants_csv(tmp_path, designed_variants_df):
    """Create a temporary designed variants CSV file."""
    csv_file = tmp_path / "mock_variants.csv"
    designed_variants_df.to_csv(csv_file, index=False)
    return csv_file


@pytest.fixture
def mock_experiment_csv(tmp_path, mock_experiment_df):
    """Create a temporary experiment CSV file."""
    csv_file = tmp_path / "mock_experiment.csv"
    mock_experiment_df.to_csv(csv_file, index=False)
    return csv_file


# =============================================================================
# Integration Test Fixtures
# =============================================================================


@pytest.fixture
def mock_config_file(tmp_path, mock_config, mock_reference_fasta, mock_variants_csv, mock_experiment_csv):
    """
    Create a complete mock configuration file with all referenced files.
    Returns path to the config YAML file.
    """
    import yaml

    # Update config paths to point to temp files
    config = mock_config.copy()
    config["ref_dir"] = str(tmp_path)
    config["reference"] = mock_reference_fasta.name
    config["variants_file"] = str(mock_variants_csv)
    config["experiment_file"] = str(mock_experiment_csv)
    config["data_dir"] = str(tmp_path)

    config_file = tmp_path / "mock_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config, f)

    return config_file
