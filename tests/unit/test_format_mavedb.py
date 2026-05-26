"""
Unit tests for workflow/rules/scripts/format_mavedb.py — the MaveDB score-set
formatter that converts a rosace or lilace score CSV into MaveDB's expected
hgvs_pro / score / sd column layout.

The invariants pinned here:

  - _BACKEND_DEFAULTS picks the right column names for each backend.
    Rosace uses `variants`/`mean`/`sd`; lilace uses `variant`/`effect`/
    `effect_se`. The two backends ship distinct column names and this
    test catches the schema drift if either backend changes.

  - normalize_hgvs converts rosace's same-AA missense notation
    (`p.(A1A)`) into MaveDB-canonical synonymous notation (`p.(A1=)`).
    mavehgvs treats these as distinct variant types, so normalization
    has to happen BEFORE mavehgvs validation.

  - Missing required columns produce a clear, actionable error that
    names the file path and lists what was available.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "workflow" / "rules" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# format_mavedb.py executes main() at import time when run via Snakemake.
# Guard against that by checking the `snakemake` global isn't present (it
# isn't in this test process), then import the helper functions directly.
import format_mavedb  # noqa: E402


# -----------------------------------------------------------------------------
# normalize_hgvs: synonymous notation conversion
# -----------------------------------------------------------------------------


class TestNormalizeHgvs:
    def test_same_aa_missense_becomes_synonymous(self):
        # Rosace emits `p.(A1A)` for synonymous variants (same-AA substitution);
        # MaveDB / mavehgvs canonical form is `p.(A1=)`. The downstream
        # mavehgvs validator treats `p.(A1A)` as a distinct (non-synonymous)
        # variant type — see the audit notes in tasks/tasks.md. Normalize
        # before validation to avoid spurious "wrong variant type" warnings
        # at MaveDB upload.
        assert format_mavedb.normalize_hgvs("p.(A1A)") == "p.(A1=)"
        assert format_mavedb.normalize_hgvs("p.(M100M)") == "p.(M100=)"

    def test_real_missense_passes_through(self):
        assert format_mavedb.normalize_hgvs("p.(A1V)") == "p.(A1V)"
        assert format_mavedb.normalize_hgvs("p.(M1L)") == "p.(M1L)"

    def test_already_canonical_synonymous_passes_through(self):
        assert format_mavedb.normalize_hgvs("p.(A1=)") == "p.(A1=)"

    def test_indel_passes_through_unchanged(self):
        # Normalize only fires on the same-AA missense shape. Anything more
        # complex (deletions, insertions, frameshifts) is delegated to
        # mavehgvs for validation in Phase 2 and not touched here.
        assert format_mavedb.normalize_hgvs("p.(A1del)") == "p.(A1del)"
        assert format_mavedb.normalize_hgvs("p.(A1_A3del)") == "p.(A1_A3del)"
        assert format_mavedb.normalize_hgvs("p.(A1_A2insGGG)") == "p.(A1_A2insGGG)"

    def test_nonsense_passes_through(self):
        assert format_mavedb.normalize_hgvs("p.(A1*)") == "p.(A1*)"


# -----------------------------------------------------------------------------
# Backend defaults: rosace and lilace ship different column names. The
# defaults must match the actual emitted schema (verified empirically in
# Phase 0 — see tasks/tasks.md audit notes).
# -----------------------------------------------------------------------------


class TestBackendDefaults:
    def test_rosace_defaults_match_emitted_schema(self):
        # Verified against results/example_experiment/rosace/cond_A_scores.csv
        # at audit time: columns are variants, ..., mean, sd, ...
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        assert d == {"hgvs_col": "variants", "score_col": "mean", "sd_col": "sd"}

    def test_lilace_defaults_match_emitted_schema(self):
        # Verified against pimentellab/lilace R/lilace.R .get_score_df():
        # singular `variant`, `effect`, `effect_se`. NOT the same as rosace.
        # If this assertion fires because someone "unified" the defaults,
        # lilace's deposit path would silently look for nonexistent columns.
        d = format_mavedb._BACKEND_DEFAULTS["lilace"]
        assert d == {
            "hgvs_col": "variant",
            "score_col": "effect",
            "sd_col": "effect_se",
        }


# -----------------------------------------------------------------------------
# End-to-end: format_scores writes the right MaveDB CSV shape for both
# backends. Use tiny in-memory fixtures so the test is hermetic.
# -----------------------------------------------------------------------------


@pytest.fixture
def rosace_scores_csv(tmp_path):
    csv = tmp_path / "rosace_scores.csv"
    pd.DataFrame(
        {
            "variants": ["p.(M1L)", "p.(M1M)", "p.(A2del)"],
            "position": [1, 1, 2],
            "wildtype": ["M", "M", "A"],
            "mutation": ["L", "M", "-"],
            "type": ["missense", "synonymous", "deletion"],
            "mean": [-0.84, 0.01, -1.20],
            "sd": [0.27, 0.05, 0.40],
            "lfsr": [0.001, 0.50, 0.005],
        }
    ).to_csv(csv, index=False)
    return csv


@pytest.fixture
def lilace_scores_csv(tmp_path):
    csv = tmp_path / "lilace_scores.csv"
    pd.DataFrame(
        {
            "variant": ["p.(M1L)", "p.(M1M)", "p.(A2del)"],
            "type": ["missense", "synonymous", "deletion"],
            "position": [1, 1, 2],
            "effect": [-0.91, 0.02, -1.18],
            "effect_se": [0.30, 0.06, 0.42],
            "lfsr": [0.002, 0.45, 0.008],
        }
    ).to_csv(csv, index=False)
    return csv


class TestFormatScores:
    def test_rosace_end_to_end(self, tmp_path, rosace_scores_csv):
        out = tmp_path / "rosace_mavedb.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(rosace_scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(out),
        )
        df = pd.read_csv(out)
        assert list(df.columns) == ["hgvs_pro", "score", "sd"]
        # Synonymous normalization fired on row 2 (M1M → M1=), then the
        # whole row was converted to MAVE-HGVS form (three-letter codes,
        # no parens). End-to-end check that both stages compose correctly.
        assert list(df["hgvs_pro"]) == ["p.Met1Leu", "p.Met1=", "p.Ala2del"]
        # Scores pass through unchanged.
        assert df["score"].tolist() == [-0.84, 0.01, -1.20]
        assert df["sd"].tolist() == [0.27, 0.05, 0.40]

    def test_lilace_end_to_end(self, tmp_path, lilace_scores_csv):
        out = tmp_path / "lilace_mavedb.csv"
        d = format_mavedb._BACKEND_DEFAULTS["lilace"]
        format_mavedb.format_scores(
            scores_path=str(lilace_scores_csv),
            backend="lilace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(out),
        )
        df = pd.read_csv(out)
        assert list(df.columns) == ["hgvs_pro", "score", "sd"]
        assert list(df["hgvs_pro"]) == ["p.Met1Leu", "p.Met1=", "p.Ala2del"]
        assert df["score"].tolist() == [-0.91, 0.02, -1.18]
        assert df["sd"].tolist() == [0.30, 0.06, 0.42]

    def test_missing_score_column_raises_named_error(self, tmp_path):
        bad = tmp_path / "missing_score.csv"
        pd.DataFrame({"variants": ["p.(A1V)"]}).to_csv(bad, index=False)
        with pytest.raises(ValueError) as excinfo:
            format_mavedb.format_scores(
                scores_path=str(bad),
                backend="rosace",
                hgvs_col="variants",
                score_col="mean",
                sd_col="sd",
                output_path=str(tmp_path / "out.csv"),
            )
        msg = str(excinfo.value)
        # Error must point at the offending file and list available cols so
        # users can either fix the CSV or override the column-name config.
        assert str(bad) in msg
        assert "mean" in msg
        assert "variants" in msg  # listed under available columns

    def test_sd_column_omitted_when_absent(self, tmp_path):
        # MaveDB's `sd` column is optional. If the score CSV doesn't have a
        # standard-deviation column (e.g. an old lilace version), the
        # formatter should write out a valid 2-column MaveDB CSV rather
        # than crashing or fabricating SD values.
        csv = tmp_path / "no_sd.csv"
        pd.DataFrame(
            {"variants": ["p.(A1V)"], "mean": [-0.5]}
        ).to_csv(csv, index=False)
        out = tmp_path / "out.csv"
        format_mavedb.format_scores(
            scores_path=str(csv),
            backend="rosace",
            hgvs_col="variants",
            score_col="mean",
            sd_col="sd",  # the column we're looking for, but it's not in the input
            output_path=str(out),
        )
        df = pd.read_csv(out)
        assert list(df.columns) == ["hgvs_pro", "score"]


# -----------------------------------------------------------------------------
# to_mave_hgvs: structural conversion from rosace/lilace's one-letter
# parenthesized HGVS into MAVE-HGVS's three-letter no-paren form.
#
# MAVE-HGVS only accepts the `p.(...)` parenthesization for the population-
# level synonymous form `p.(=)`. For everything else, parens are rejected
# and the AA codes must be three-letter. `Ter` is the stop code, NOT `*`.
# These tests pin the conversion table for every variant class dumpling
# emits.
# -----------------------------------------------------------------------------


class TestToMaveHgvs:
    def test_missense(self):
        # p.(M1L) -> p.Met1Leu: outer parens stripped, one-letter expanded.
        assert format_mavedb.to_mave_hgvs("p.(M1L)") == "p.Met1Leu"
        assert format_mavedb.to_mave_hgvs("p.(A100V)") == "p.Ala100Val"

    def test_synonymous_at_known_position(self):
        # p.(M1=) -> p.Met1=: the WT one-letter becomes three-letter, the
        # `=` stays. Position MUST be preserved (it's the load-bearing
        # difference from the population-level p.(=) form).
        assert format_mavedb.to_mave_hgvs("p.(M1=)") == "p.Met1="
        assert format_mavedb.to_mave_hgvs("p.(C22=)") == "p.Cys22="

    def test_population_level_synonymous_passes_through(self):
        # p.(=) is the one parenthesized form MAVE-HGVS accepts. Don't
        # strip its parens. Dumpling doesn't produce this today but the
        # contract should be explicit so future code can.
        assert format_mavedb.to_mave_hgvs("p.(=)") == "p.(=)"

    def test_single_residue_deletion(self):
        assert format_mavedb.to_mave_hgvs("p.(A1del)") == "p.Ala1del"
        assert format_mavedb.to_mave_hgvs("p.(G18del)") == "p.Gly18del"

    def test_multi_residue_deletion(self):
        # Both endpoints carry their WT amino acid. Standard HGVS form.
        # MAVE-HGVS expects `p.<wt3>{start}_<wt3>{end}del`.
        assert (
            format_mavedb.to_mave_hgvs("p.(A1_C3del)") == "p.Ala1_Cys3del"
        )
        assert (
            format_mavedb.to_mave_hgvs("p.(Q7_N19del)") == "p.Gln7_Asn19del"
        )

    def test_insertion(self):
        # The inserted residues are themselves one-letter coded in the
        # rosace output; expand each. Endpoints get the same one→three
        # treatment as deletions.
        assert (
            format_mavedb.to_mave_hgvs("p.(A1_A2insGGG)")
            == "p.Ala1_Ala2insGlyGlyGly"
        )
        assert (
            format_mavedb.to_mave_hgvs("p.(H7_Q8insS)")
            == "p.His7_Gln8insSer"
        )

    def test_nonsense(self):
        # `*` must become `Ter`; MAVE-HGVS explicitly rejects `*` as the
        # stop code (audit finding, see tasks/tasks.md).
        assert format_mavedb.to_mave_hgvs("p.(A1*)") == "p.Ala1Ter"
        assert format_mavedb.to_mave_hgvs("p.(M1*)") == "p.Met1Ter"

    def test_unknown_shape_raises_actionable(self):
        # The error message must list every supported shape so users (or
        # future me) reading the traceback know what dumpling SHOULD have
        # produced. If a new HGVS shape lands in scoring output, the
        # converter needs to grow a new pattern — this test forces the
        # decision to be visible rather than silent.
        with pytest.raises(ValueError) as excinfo:
            format_mavedb.to_mave_hgvs("p.(M1?)")
        msg = str(excinfo.value)
        assert "Unrecognized HGVS shape" in msg
        assert "missense" in msg
        assert "insertion" in msg


# -----------------------------------------------------------------------------
# End-to-end format_scores with MAVE-HGVS conversion applied. Rosace's
# one-letter parens go in, three-letter no-parens come out.
# -----------------------------------------------------------------------------


class TestFormatScoresMaveHgvs:
    def test_rosace_output_is_mave_hgvs(self, tmp_path):
        csv = tmp_path / "rosace_scores.csv"
        pd.DataFrame(
            {
                "variants": [
                    "p.(M1L)",  # missense
                    "p.(M1M)",  # rosace same-AA → should normalize then convert
                    "p.(A1del)",  # single deletion
                    "p.(A1_C3del)",  # multi deletion
                    "p.(A1*)",  # nonsense
                ],
                "mean": [-0.84, 0.01, -1.20, -1.55, -2.10],
                "sd": [0.27, 0.05, 0.40, 0.45, 0.60],
            }
        ).to_csv(csv, index=False)
        out = tmp_path / "rosace_mavedb.csv"

        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(out),
        )
        df = pd.read_csv(out)
        # No parens on substitutions/deletions/nonsense; three-letter
        # everywhere; Ter for stop; synonymous keeps its position.
        assert list(df["hgvs_pro"]) == [
            "p.Met1Leu",
            "p.Met1=",
            "p.Ala1del",
            "p.Ala1_Cys3del",
            "p.Ala1Ter",
        ]


# -----------------------------------------------------------------------------
# mavehgvs validation hook. The library may not be installed in the dev
# environment (it's pip-only and not always available), so use importorskip
# to gate the live-validation tests. The hook's no-op-when-missing behavior
# is tested separately without needing the library.
# -----------------------------------------------------------------------------


class TestMaveHgvsValidation:
    def test_no_op_when_library_missing(self, monkeypatch):
        # If mavehgvs isn't importable, the validator should pass strings
        # through unchanged rather than crashing. Simulate by patching out
        # the import-time global.
        monkeypatch.setattr(format_mavedb, "_MaveVariant", None)
        assert (
            format_mavedb._validate_mave_hgvs("p.Met1Leu", 0) == "p.Met1Leu"
        )
        # Even an obviously-bad string passes through when validation
        # can't run — the no-op contract is unconditional, not "safe
        # strings only".
        assert (
            format_mavedb._validate_mave_hgvs("garbage", 0) == "garbage"
        )

    def test_validation_fires_when_library_present(self):
        # Run the real validator path. If mavehgvs is installed in the env,
        # a known-good MAVE-HGVS string round-trips; otherwise skip.
        pytest.importorskip(
            "mavehgvs",
            reason="mavehgvs not installed in this env; CI exercises this path",
        )
        # Reload the symbol in case it was monkey-patched away by a sibling
        # test in this session.
        from mavehgvs import Variant as MaveVariant

        format_mavedb._MaveVariant = MaveVariant
        assert (
            format_mavedb._validate_mave_hgvs("p.Met1Leu", 0) == "p.Met1Leu"
        )

    def test_validation_error_names_row_and_value(self):
        pytest.importorskip(
            "mavehgvs",
            reason="mavehgvs not installed in this env; CI exercises this path",
        )
        from mavehgvs import Variant as MaveVariant

        format_mavedb._MaveVariant = MaveVariant
        # `p.(M1L)` (one-letter parenthesized) is exactly what mavehgvs
        # rejects — that's the failure mode this whole refactor exists to
        # surface. The error message must point at the row index and the
        # offending value so a user with a 10k-row CSV can find it.
        with pytest.raises(ValueError) as excinfo:
            format_mavedb._validate_mave_hgvs("p.(M1L)", 42)
        msg = str(excinfo.value)
        assert "Row 42" in msg
        assert "p.(M1L)" in msg
