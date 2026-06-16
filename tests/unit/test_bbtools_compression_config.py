"""
Schema- and shim-level tests for the `bbtools_compression` configuration knob,
which replaced the deprecated `bbtools_use_bgzip` bool to enable parallel pigz.

The invariants pinned here:

  - bbtools_compression is not in `required`, so old configs (pre-pigz) without
    the key still validate.
  - the default is "pigz", so users who specified neither key get the faster
    parallelized (de)compression path automatically.
  - the enum is {pigz, bgzip, none}, so typos and accidental other values
    are rejected at schema-validation time.
  - snakemake.utils.validate injects the default at runtime — which is
    the mechanism common.smk actually relies on for backward compat with
    configs that mention neither key.
  - the translate_legacy_bbtools_compression shim preserves byte-compatible
    behavior for old configs:
      bbtools_use_bgzip: true  →  bbtools_compression: bgzip
      bbtools_use_bgzip: false →  bbtools_compression: none

Same gating pattern as test_aligner_config.py: skip cleanly if jsonschema
or snakemake.utils isn't importable.
"""

import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "workflow" / "schemas" / "config.schema.yaml"
SCRIPTS_DIR = REPO_ROOT / "workflow" / "rules" / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from script_utils import translate_legacy_bbtools_compression  # noqa: E402


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
# Schema shape: the document on disk has the expected `bbtools_compression`
# block. Runs with no extra deps and catches the most common drift (someone
# renames the property, drops the enum, etc.).
# -----------------------------------------------------------------------------


class TestSchemaShape:
    def test_bbtools_compression_property_present(self):
        schema = _load_schema()
        assert "bbtools_compression" in schema["properties"]
        prop = schema["properties"]["bbtools_compression"]
        assert prop["type"] == "string"
        assert set(prop["enum"]) == {"pigz", "bgzip", "none"}
        assert prop["default"] == "pigz"

    def test_bbtools_compression_not_required(self):
        """Old configs from before pigz lands don't have this key. The schema
        must not require it — otherwise users on the old config shape would be
        locked out by validation."""
        schema = _load_schema()
        assert "bbtools_compression" not in schema.get("required", [])

    def test_legacy_bbtools_use_bgzip_not_in_schema(self):
        """The deprecated key is intentionally absent from the schema.
        jsonschema's default `additionalProperties: true` lets old configs
        through, and the runtime shim pops the key before downstream code
        reads it. If someone ever re-adds it to the schema, the deprecation
        path becomes ambiguous (two valid keys in the validated config)."""
        schema = _load_schema()
        assert "bbtools_use_bgzip" not in schema["properties"]


# -----------------------------------------------------------------------------
# Plain jsonschema validation: enum enforcement, type checks. Doesn't fill
# defaults (jsonschema doesn't by default).
# -----------------------------------------------------------------------------


jsonschema = pytest.importorskip(
    "jsonschema",
    reason="jsonschema not installed; skipping schema enum tests",
)


class TestEnumEnforcement:
    def test_bbtools_compression_pigz_accepted(self):
        config = _minimal_valid_config()
        config["bbtools_compression"] = "pigz"
        jsonschema.validate(config, _load_schema())

    def test_bbtools_compression_bgzip_accepted(self):
        config = _minimal_valid_config()
        config["bbtools_compression"] = "bgzip"
        jsonschema.validate(config, _load_schema())

    def test_bbtools_compression_none_accepted(self):
        config = _minimal_valid_config()
        config["bbtools_compression"] = "none"
        jsonschema.validate(config, _load_schema())

    def test_bbtools_compression_unknown_rejected(self):
        config = _minimal_valid_config()
        config["bbtools_compression"] = "zstd"
        with pytest.raises(jsonschema.ValidationError) as excinfo:
            jsonschema.validate(config, _load_schema())
        msg = str(excinfo.value)
        assert "zstd" in msg or "enum" in msg

    def test_bbtools_compression_omitted_validates(self):
        """Pre-pigz configs don't set this key. Schema must accept that shape;
        defaulting happens at runtime via snakemake.utils.validate."""
        config = _minimal_valid_config()
        assert "bbtools_compression" not in config
        jsonschema.validate(config, _load_schema())

    def test_legacy_bbtools_use_bgzip_validates(self):
        """The deprecated key is not in the schema; jsonschema's default
        permissive policy on extra properties means a legacy config still
        passes validation (the runtime shim handles translation). If the
        schema ever flips to additionalProperties=false, the shim must run
        before validation."""
        config = _minimal_valid_config()
        config["bbtools_use_bgzip"] = True
        jsonschema.validate(config, _load_schema())


# -----------------------------------------------------------------------------
# Default injection via snakemake.utils.validate. This is the runtime behavior
# common.smk actually relies on, so locking it in matters for backward compat:
# if a user's existing config doesn't have `bbtools_compression`, they MUST
# end up on the pigz path (the new default) automatically.
# -----------------------------------------------------------------------------


class TestDefaultInjection:
    @pytest.fixture(autouse=True)
    def _require_snakemake(self):
        self.snakemake_validate = pytest.importorskip(
            "snakemake.utils",
            reason="snakemake not installed; skipping default-injection tests",
        ).validate

    def test_bbtools_compression_default_is_pigz(self):
        """A config that mentions neither key (the modern default-path case)
        must end up on pigz after validation. If this assertion fires, the
        new default has drifted."""
        config = _minimal_valid_config()
        assert "bbtools_compression" not in config
        assert "bbtools_use_bgzip" not in config
        self.snakemake_validate(config, str(SCHEMA_PATH))
        assert config.get("bbtools_compression") == "pigz"


# -----------------------------------------------------------------------------
# Back-compat shim: translate_legacy_bbtools_compression mutates the config
# in place and returns warning strings. Byte-compatibility with the prior
# helper is the load-bearing invariant — old configs must keep producing the
# same BBTools flag string after translation.
# -----------------------------------------------------------------------------


class TestLegacyShim:
    def test_legacy_true_translates_to_bgzip(self):
        """Old default `bbtools_use_bgzip: true` produced an empty flag
        string, which is what the new `bgzip` mode does. Translating to
        `bgzip` preserves byte-identical BBTools invocations."""
        config = {"bbtools_use_bgzip": True}
        warnings = translate_legacy_bbtools_compression(config)
        assert config.get("bbtools_compression") == "bgzip"
        assert "bbtools_use_bgzip" not in config
        assert len(warnings) == 1
        assert "deprecated" in warnings[0]

    def test_legacy_false_translates_to_none(self):
        """Old `bbtools_use_bgzip: false` produced `bgzip=f unbgzip=f`,
        which is what the new `none` mode does."""
        config = {"bbtools_use_bgzip": False}
        warnings = translate_legacy_bbtools_compression(config)
        assert config.get("bbtools_compression") == "none"
        assert "bbtools_use_bgzip" not in config
        assert len(warnings) == 1
        assert "deprecated" in warnings[0]

    def test_absent_legacy_key_no_op(self):
        """If the legacy key isn't present, the shim is a no-op and emits
        no warnings — the modern default path."""
        config = {}
        warnings = translate_legacy_bbtools_compression(config)
        assert warnings == []
        assert "bbtools_compression" not in config

    def test_both_keys_set_new_wins(self):
        """If a user has both keys set (unusual, but possible after a
        copy-paste from old docs), the new key wins and we warn about the
        redundant legacy one. Silently dropping the legacy choice would be
        a worse outcome than the warning."""
        config = {"bbtools_use_bgzip": False, "bbtools_compression": "pigz"}
        warnings = translate_legacy_bbtools_compression(config)
        assert config["bbtools_compression"] == "pigz"
        assert "bbtools_use_bgzip" not in config
        assert len(warnings) == 1
        assert "ignoring" in warnings[0]
