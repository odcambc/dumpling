"""
Schema-level tests for the `lilace_seed` configuration knob (Determinism
Phase 1). The invariants pinned here:

  - lilace_seed is not in `required`, so old configs (pre-seed) without
    the key still validate.
  - the default is null (None), which means "auto-generate at
    rule-invocation time and log" — preserves prior behavior for users
    who upgrade without touching their config.
  - both null and positive integers are accepted; arbitrary strings are
    rejected.
  - snakemake.utils.validate injects null at runtime, which the R script
    detects via is.null() to decide between configured and auto seeds.
  - there is intentionally no `rosace_seed` knob. RunRosace doesn't
    accept a seed argument (see pimentellab/rosace runROSACE.R); adding
    a knob would silently no-op and mislead users into thinking they
    were getting reproducible rosace chains when seed=100 is already
    hardcoded.

Same gating pattern as test_aligner_config.py and
test_bbtools_compression_config.py.
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


class TestSchemaShape:
    def test_lilace_seed_property_present(self):
        schema = _load_schema()
        assert "lilace_seed" in schema["properties"]
        prop = schema["properties"]["lilace_seed"]
        # Type is the YAML list [integer, "null"]
        assert set(prop["type"]) == {"integer", "null"}
        assert prop["default"] is None

    def test_lilace_seed_not_required(self):
        schema = _load_schema()
        assert "lilace_seed" not in schema.get("required", [])

    def test_no_rosace_seed_knob(self):
        """RunRosace doesn't accept a seed argument upstream; a config
        knob would silently no-op. The absence of `rosace_seed` from
        the schema is intentional — see pimentellab/rosace
        R/runROSACE.R (MCMCRunStan hardcodes seed=100)."""
        schema = _load_schema()
        assert "rosace_seed" not in schema["properties"]


jsonschema = pytest.importorskip(
    "jsonschema",
    reason="jsonschema not installed; skipping schema enum tests",
)


class TestEnumEnforcement:
    def test_lilace_seed_null_accepted(self):
        config = _minimal_valid_config()
        config["lilace_seed"] = None
        jsonschema.validate(config, _load_schema())

    def test_lilace_seed_integer_accepted(self):
        config = _minimal_valid_config()
        config["lilace_seed"] = 42
        jsonschema.validate(config, _load_schema())

    def test_lilace_seed_zero_accepted(self):
        """Stan accepts 0 as a valid seed; the schema shouldn't reject
        it. Users who want a deterministic 'no-salt' run will reach
        for 0."""
        config = _minimal_valid_config()
        config["lilace_seed"] = 0
        jsonschema.validate(config, _load_schema())

    def test_lilace_seed_string_rejected(self):
        config = _minimal_valid_config()
        config["lilace_seed"] = "deadbeef"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(config, _load_schema())

    def test_lilace_seed_omitted_validates(self):
        config = _minimal_valid_config()
        assert "lilace_seed" not in config
        jsonschema.validate(config, _load_schema())


class TestDefaultInjection:
    @pytest.fixture(autouse=True)
    def _require_snakemake(self):
        self.snakemake_validate = pytest.importorskip(
            "snakemake.utils",
            reason="snakemake not installed; skipping default-injection tests",
        ).validate

    def test_lilace_seed_default_is_null(self):
        """Default lilace_seed must be null — that's what
        run_lilace.R's resolve_lilace_seed checks with is.null() to
        decide between auto-generate vs configured."""
        config = _minimal_valid_config()
        assert "lilace_seed" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("lilace_seed") is None
