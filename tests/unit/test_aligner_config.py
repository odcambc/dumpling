"""
Schema-level tests for the `aligner` configuration knob introduced when
minimap2 was added as an opt-in alternative to bbmap.

The invariants pinned here:

  - aligner is not in `required`, so old configs (pre-minimap2) without
    the key still validate.
  - the default is "bbmap", so old configs run unchanged on the bbmap path.
  - the enum is {bbmap, minimap2}, so typos and accidental other values
    are rejected at schema-validation time.
  - snakemake.utils.validate injects the default at runtime — which is
    the mechanism `common.smk` actually relies on for backward compat.

Same gating pattern as test_scoring_backend_config.py: skip cleanly if
jsonschema or snakemake.utils isn't importable.
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
# Schema shape: the document on disk has the expected `aligner` block.
# Runs with no extra deps and catches the most common drift (someone renames
# the property, drops the enum, etc.).
# -----------------------------------------------------------------------------


class TestSchemaShape:
    def test_aligner_property_present(self):
        schema = _load_schema()
        assert "aligner" in schema["properties"]
        prop = schema["properties"]["aligner"]
        assert prop["type"] == "string"
        assert set(prop["enum"]) == {"bbmap", "minimap2"}
        assert prop["default"] == "bbmap"

    def test_aligner_not_required(self):
        """Old configs from before minimap2 lands don't have an `aligner`
        key. The schema must not require it — otherwise users on the old
        config shape would be locked out by validation."""
        schema = _load_schema()
        assert "aligner" not in schema.get("required", [])


# -----------------------------------------------------------------------------
# Plain jsonschema validation: enum enforcement, type checks. Doesn't fill
# defaults (jsonschema doesn't by default).
# -----------------------------------------------------------------------------


jsonschema = pytest.importorskip(
    "jsonschema",
    reason="jsonschema not installed; skipping schema enum tests",
)


class TestEnumEnforcement:
    def test_aligner_bbmap_accepted(self):
        config = _minimal_valid_config()
        config["aligner"] = "bbmap"
        jsonschema.validate(config, _load_schema())

    def test_aligner_minimap2_accepted(self):
        config = _minimal_valid_config()
        config["aligner"] = "minimap2"
        jsonschema.validate(config, _load_schema())

    def test_aligner_unknown_rejected(self):
        config = _minimal_valid_config()
        config["aligner"] = "bowtie2"
        with pytest.raises(jsonschema.ValidationError) as excinfo:
            jsonschema.validate(config, _load_schema())
        msg = str(excinfo.value)
        assert "bowtie2" in msg or "enum" in msg

    def test_aligner_omitted_validates(self):
        """Pre-minimap2 configs don't set this key. Schema must accept that
        shape — defaulting happens at runtime in common.smk via setdefault."""
        config = _minimal_valid_config()
        assert "aligner" not in config
        jsonschema.validate(config, _load_schema())


# -----------------------------------------------------------------------------
# Default injection via snakemake.utils.validate. This is the runtime behavior
# common.smk actually relies on, so locking it in matters for backward compat:
# if a user's existing config doesn't have `aligner`, they MUST end up on the
# bbmap path with no behavior change.
# -----------------------------------------------------------------------------


class TestDefaultInjection:
    @pytest.fixture(autouse=True)
    def _require_snakemake(self):
        self.snakemake_validate = pytest.importorskip(
            "snakemake.utils",
            reason="snakemake not installed; skipping default-injection tests",
        ).validate

    def test_aligner_default_is_bbmap(self):
        """The headline backward-compat invariant: a config that doesn't
        mention `aligner` must run on the bbmap path. If this assertion
        fires, an existing user's pipeline is about to switch aligners
        out from under them."""
        config = _minimal_valid_config()
        assert "aligner" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("aligner") == "bbmap"
