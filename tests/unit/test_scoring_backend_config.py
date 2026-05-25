"""
Schema-level tests for the scoring_backend / lilace_local / mem_lilace
configuration knobs introduced for the Lilace scoring backend.

These tests cover the config.schema.yaml contract and the
snakemake-utils-validate default-injection behavior. They are gated on
the availability of snakemake (for default injection) or jsonschema (for
plain validation), and skip cleanly when neither is installed.
"""

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "workflow" / "schemas" / "config.schema.yaml"


def _load_schema():
    with open(SCHEMA_PATH) as fh:
        return yaml.safe_load(fh)


def _minimal_valid_config():
    """A config that satisfies all `required` fields in config.schema.yaml.

    Intentionally omits all optional fields so the tests can observe whether
    defaults are filled or whether enums are enforced.
    """
    return {
        "data_dir": "example/data",
        "experiment": "test_experiment",
        "experiment_file": "config/example.csv",
        "ref_dir": "references",
        "reference": "example_ref.fasta",
        "orf": "1-300",
        "adapters": "resources/adapters.fa",
        "contaminants": ["resources/sequencing_artifacts.fa.gz"],
    }


# -----------------------------------------------------------------------------
# Schema shape: enums and properties are present in the document on disk.
# These run with no extra deps and catch the most common drift.
# -----------------------------------------------------------------------------


class TestSchemaShape:
    def test_scoring_backend_property_present(self):
        schema = _load_schema()
        assert "scoring_backend" in schema["properties"]
        prop = schema["properties"]["scoring_backend"]
        assert prop["type"] == "string"
        assert set(prop["enum"]) == {"rosace", "lilace"}
        assert prop["default"] == "rosace"

    def test_mem_lilace_property_present(self):
        schema = _load_schema()
        assert "mem_lilace" in schema["properties"]
        prop = schema["properties"]["mem_lilace"]
        assert prop["type"] == "integer"
        assert prop["default"] == 16000

    def test_lilace_local_property_present(self):
        schema = _load_schema()
        assert "lilace_local" in schema["properties"]
        prop = schema["properties"]["lilace_local"]
        assert prop["type"] == "boolean"
        assert prop["default"] is False


# -----------------------------------------------------------------------------
# Plain jsonschema validation: enum enforcement, type checks. No default
# injection (jsonschema doesn't fill defaults by default).
# -----------------------------------------------------------------------------


jsonschema = pytest.importorskip(
    "jsonschema",
    reason="jsonschema not installed; skipping schema enum tests",
)


class TestEnumEnforcement:
    def test_scoring_backend_rosace_accepted(self):
        config = _minimal_valid_config()
        config["scoring_backend"] = "rosace"
        jsonschema.validate(config, _load_schema())

    def test_scoring_backend_lilace_accepted(self):
        config = _minimal_valid_config()
        config["scoring_backend"] = "lilace"
        jsonschema.validate(config, _load_schema())

    def test_scoring_backend_unknown_rejected(self):
        config = _minimal_valid_config()
        config["scoring_backend"] = "enrich2"
        with pytest.raises(jsonschema.ValidationError) as excinfo:
            jsonschema.validate(config, _load_schema())
        # Error message should reference the offending value or the enum.
        msg = str(excinfo.value)
        assert "enrich2" in msg or "enum" in msg

    def test_scoring_backend_omitted_validates(self):
        """Omitting scoring_backend is fine — it has a default and is not required."""
        config = _minimal_valid_config()
        config.pop("scoring_backend", None)
        jsonschema.validate(config, _load_schema())


# -----------------------------------------------------------------------------
# Default injection via snakemake.utils.validate. This is the behavior
# common.smk actually relies on, so worth covering even though it requires
# snakemake to be importable.
# -----------------------------------------------------------------------------


class TestDefaultInjection:
    """Verify snakemake.utils.validate fills in schema defaults."""

    @pytest.fixture(autouse=True)
    def _require_snakemake(self):
        self.snakemake_validate = pytest.importorskip(
            "snakemake.utils",
            reason="snakemake not installed; skipping default-injection tests",
        ).validate

    def test_scoring_backend_default_is_rosace(self):
        config = _minimal_valid_config()
        assert "scoring_backend" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("scoring_backend") == "rosace"

    def test_mem_lilace_default_is_16000(self):
        config = _minimal_valid_config()
        assert "mem_lilace" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("mem_lilace") == 16000

    def test_lilace_local_default_is_false(self):
        config = _minimal_valid_config()
        assert "lilace_local" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("lilace_local") is False
