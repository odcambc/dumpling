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
        # MAVE-HGVS explicitly rejects `*` as the stop code; canonical
        # is `Ter`. Dumpling internally uses `X` for stops (per
        # process_variants.aa_3to1_dict's "Stp" -> "X" mapping) so that
        # is what shows up in real rosace/lilace score CSV `variants`
        # columns -- `p.(A1X)`, not `p.(A1*)`. Empirically verified
        # against results/example_experiment/rosace/cond_A_scores.csv
        # 2026-05-27: 78 rows of the form `p.([A-Z][0-9]+X)`, all
        # convert cleanly to `p.<wt3>{pos}Ter`. We still accept `*` as
        # an alternate notation for robustness.
        assert format_mavedb.to_mave_hgvs("p.(A1X)") == "p.Ala1Ter"
        assert format_mavedb.to_mave_hgvs("p.(M1X)") == "p.Met1Ter"
        assert format_mavedb.to_mave_hgvs("p.(S2X)") == "p.Ser2Ter"
        # Alternate stop notation (kept as fallback path).
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


# -----------------------------------------------------------------------------
# main()'s snakemake-vs-CLI dispatch. The original guard used
# `"snakemake" in dir()` inside main(), but dir() returns LOCAL names of
# the function — not module globals — so the snakemake-injected global
# was never visible, and main() silently fell through to the argparse path.
# Surfaced empirically on 2026-05-27 when wiring deposit_to_mavedb default-on
# into get_input and running the rule via snakemake: the script crashed
# with "argparse: error: the following arguments are required:
# --scores, --output" because sys.argv didn't contain those.
#
# Lock in the correct dispatch via globals().
# -----------------------------------------------------------------------------


class TestMainSnakemakeDispatch:
    def test_main_takes_snakemake_path_when_snakemake_global_present(
        self, monkeypatch, tmp_path
    ):
        from types import SimpleNamespace

        csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {"variants": ["p.(A1V)", "p.(M1X)"], "mean": [-0.5, -1.2], "sd": [0.1, 0.2]}
        ).to_csv(csv, index=False)
        out = tmp_path / "out.csv"

        # Inject a snakemake-shaped object into the module's globals,
        # mirroring what Snakemake's `script:` directive does at runtime.
        # If main() looked at dir() instead of globals() (the bug we just
        # fixed), this attribute would be invisible to it and main()
        # would fall through to argparse + crash on missing args.
        # snakemake.output is now a NamedTuple-like namespace because the
        # rule declares two outputs (scores + counts). When counts isn't
        # declared (this test exercises the scores-only path), only .scores
        # is set on the namespace.
        mock_snakemake = SimpleNamespace(
            config={"mavedb": {}},
            params=SimpleNamespace(backend="rosace"),
            input=SimpleNamespace(scores=str(csv)),
            output=SimpleNamespace(scores=str(out)),
        )
        monkeypatch.setattr(
            format_mavedb, "snakemake", mock_snakemake, raising=False
        )

        format_mavedb.main()

        # If we got here, main() took the snakemake path. Verify the
        # output landed where snakemake.output[0] pointed.
        assert out.exists()
        df = pd.read_csv(out)
        assert list(df.columns) == ["hgvs_pro", "score", "sd"]
        # X stop code → Ter (regression-coverage doubled with the
        # earlier TestToMaveHgvs.test_nonsense).
        assert "p.Ala1Val" in df["hgvs_pro"].tolist()
        assert "p.Met1Ter" in df["hgvs_pro"].tolist()

    def test_main_falls_through_to_argparse_when_no_snakemake_global(
        self, monkeypatch
    ):
        # Belt-and-braces: when no `snakemake` global is present, main()
        # MUST go down the CLI path (which argparse-errors out without
        # required args). Catches a regression where someone "fixes" the
        # guard by removing the snakemake-mode detection entirely.
        monkeypatch.setattr(sys, "argv", ["format_mavedb.py"])
        # Ensure no leftover snakemake attribute from a sibling test.
        if hasattr(format_mavedb, "snakemake"):
            monkeypatch.delattr(format_mavedb, "snakemake")
        with pytest.raises(SystemExit):
            format_mavedb.main()


# -----------------------------------------------------------------------------
# _aggregate_sample_counts: read one processed_counts/{sample}.csv and sum
# across the codons that encode the same protein variant. dumpling emits one
# row per (variant, codon), so M1F from TTC and M1F from TTT must collapse
# to a single MAVE-HGVS key with the combined count.
# -----------------------------------------------------------------------------


class TestAggregateSampleCounts:
    def test_multi_codon_aggregates_to_protein_variant(self, tmp_path):
        # Same protein variant from two codons must sum into one row of the
        # output series. If this regresses, MaveDB count tables will
        # double-count protein variants whose codons happen to land on the
        # same key after MAVE-HGVS conversion.
        csv = tmp_path / "A_R1_T0.csv"
        pd.DataFrame(
            {
                "count": [10, 5, 200],
                "hgvs": ["p.(M1F)", "p.(M1F)", "p.(M1I)"],
            }
        ).to_csv(csv, index=False)
        result = format_mavedb._aggregate_sample_counts(str(csv))
        assert result["p.Met1Phe"] == 15
        assert result["p.Met1Ile"] == 200
        assert len(result) == 2

    def test_synonymous_codons_collapse_to_canonical_key(self, tmp_path):
        # rosace's per-codon hgvs for a synonymous variant is "p.(X1X)"
        # (codon-specific same-AA notation). normalize_hgvs rewrites these
        # to canonical "p.(X1=)" before MAVE-HGVS conversion, so multiple
        # codons encoding the same WT at the same position collapse to a
        # single "p.<wt3>1=" key.
        csv = tmp_path / "sample.csv"
        pd.DataFrame(
            {
                "count": [50, 30, 20],
                "hgvs": ["p.(M1M)", "p.(M1M)", "p.(M1M)"],
            }
        ).to_csv(csv, index=False)
        result = format_mavedb._aggregate_sample_counts(str(csv))
        assert result["p.Met1="] == 100
        assert len(result) == 1

    def test_missing_required_columns_raises_named_error(self, tmp_path):
        # Error must name the file path AND the missing column so a user
        # with many per-sample CSVs can find the bad one quickly.
        csv = tmp_path / "broken.csv"
        pd.DataFrame({"hgvs": ["p.(A1V)"]}).to_csv(csv, index=False)
        with pytest.raises(ValueError) as excinfo:
            format_mavedb._aggregate_sample_counts(str(csv))
        msg = str(excinfo.value)
        assert str(csv) in msg
        assert "count" in msg


# -----------------------------------------------------------------------------
# format_scores with counts emission: union variant set across scores +
# counts, NaN-pad scores for variants only in counts, zero-pad counts for
# variants only in scores, and preserve sample-name column order. The
# MaveDB spec requires score and count tables to share the same variant set;
# violating that would cause deposit-time validation errors.
# -----------------------------------------------------------------------------


@pytest.fixture
def two_sample_counts(tmp_path):
    """Two per-sample count CSVs with overlapping but not identical variant
    sets. A_R1_T0 has M1F + M1I + M1=; A_R1_T1 has M1F + M1I + M1= + A2V."""
    sample_a = tmp_path / "A_R1_T0.csv"
    pd.DataFrame(
        {
            "count": [10, 5, 200, 50],
            "hgvs": ["p.(M1F)", "p.(M1F)", "p.(M1I)", "p.(M1M)"],
        }
    ).to_csv(sample_a, index=False)

    sample_b = tmp_path / "A_R1_T1.csv"
    pd.DataFrame(
        {
            "count": [12, 180, 40, 999],
            "hgvs": ["p.(M1F)", "p.(M1I)", "p.(M1M)", "p.(A2V)"],
        }
    ).to_csv(sample_b, index=False)

    return [sample_a, sample_b], ["A_R1_T0", "A_R1_T1"]


class TestFormatScoresWithCounts:
    def test_score_and_count_csvs_share_variant_set(
        self, tmp_path, two_sample_counts
    ):
        counts_paths, sample_names = two_sample_counts
        scores_csv = tmp_path / "scores.csv"
        # Score CSV has M1F + M1= only (rosace filtered M1I and A2V out).
        pd.DataFrame(
            {
                "variants": ["p.(M1F)", "p.(M1M)"],
                "mean": [-0.5, 0.02],
                "sd": [0.1, 0.05],
            }
        ).to_csv(scores_csv, index=False)

        score_out = tmp_path / "out_score.csv"
        count_out = tmp_path / "out_count.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(score_out),
            counts_paths=[str(p) for p in counts_paths],
            sample_names=sample_names,
            counts_output_path=str(count_out),
        )
        score_df = pd.read_csv(score_out)
        count_df = pd.read_csv(count_out)

        # MaveDB spec requirement: same variants, same row order. Without
        # this, the deposit-time validator rejects the upload.
        assert score_df["hgvs_pro"].tolist() == count_df["hgvs_pro"].tolist()
        assert set(score_df["hgvs_pro"]) == {
            "p.Met1Phe",
            "p.Met1=",
            "p.Met1Ile",
            "p.Ala2Val",
        }

    def test_variant_only_in_counts_gets_nan_score(
        self, tmp_path, two_sample_counts
    ):
        # A2V is only in counts (filtered out by rosace). The score CSV
        # must carry that row with NaN score/sd so depositors see the
        # filter outcomes explicitly rather than silently losing signal.
        counts_paths, sample_names = two_sample_counts
        scores_csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {
                "variants": ["p.(M1F)", "p.(M1M)"],
                "mean": [-0.5, 0.02],
                "sd": [0.1, 0.05],
            }
        ).to_csv(scores_csv, index=False)
        score_out = tmp_path / "out_score.csv"
        count_out = tmp_path / "out_count.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(score_out),
            counts_paths=[str(p) for p in counts_paths],
            sample_names=sample_names,
            counts_output_path=str(count_out),
        )
        score_df = pd.read_csv(score_out).set_index("hgvs_pro")
        # Scored variants keep their values.
        assert score_df.loc["p.Met1Phe", "score"] == -0.5
        # Counts-only variants are NaN-padded.
        assert pd.isna(score_df.loc["p.Met1Ile", "score"])
        assert pd.isna(score_df.loc["p.Met1Ile", "sd"])
        assert pd.isna(score_df.loc["p.Ala2Val", "score"])

    def test_variant_only_in_scores_gets_zero_count(self, tmp_path):
        # Symmetric edge case: a variant scored but absent from a sample's
        # count file gets 0 (not NaN) in that sample's count column. In
        # practice scoring takes counts as input so this is rare, but the
        # union path treats both sides the same way and the behavior
        # should be explicit.
        sample = tmp_path / "S.csv"
        pd.DataFrame({"count": [10], "hgvs": ["p.(M1F)"]}).to_csv(
            sample, index=False
        )
        scores_csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {
                "variants": ["p.(M1F)", "p.(M1L)"],
                "mean": [-0.5, -0.8],
                "sd": [0.1, 0.2],
            }
        ).to_csv(scores_csv, index=False)

        score_out = tmp_path / "out_score.csv"
        count_out = tmp_path / "out_count.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(score_out),
            counts_paths=[str(sample)],
            sample_names=["S"],
            counts_output_path=str(count_out),
        )
        count_df = pd.read_csv(count_out).set_index("hgvs_pro")
        assert count_df.loc["p.Met1Phe", "S"] == 10
        # Scored-only variant is zero, not NaN — fillna(0).astype(int) at
        # construction time guarantees the count column stays integer.
        assert count_df.loc["p.Met1Leu", "S"] == 0

    def test_sample_column_order_matches_input(self, tmp_path):
        # The Snakemake rule passes a deliberate sample order (sorted by
        # replicate, time). Output column order MUST match the caller's
        # order — depositors interpret column position as time-course
        # progression. Use a non-alphabetical order to catch any silent
        # alphabetization.
        order = ["A_R2_T0", "A_R1_T0", "A_R1_T1"]
        paths = []
        for s in order:
            p = tmp_path / f"{s}.csv"
            pd.DataFrame({"count": [1], "hgvs": ["p.(M1F)"]}).to_csv(
                p, index=False
            )
            paths.append(p)

        scores_csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {"variants": ["p.(M1F)"], "mean": [-0.5], "sd": [0.1]}
        ).to_csv(scores_csv, index=False)

        score_out = tmp_path / "out_score.csv"
        count_out = tmp_path / "out_count.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(score_out),
            counts_paths=[str(p) for p in paths],
            sample_names=order,
            counts_output_path=str(count_out),
        )
        count_df = pd.read_csv(count_out)
        assert list(count_df.columns) == ["hgvs_pro"] + order

    def test_misaligned_sample_names_raises(self, tmp_path, two_sample_counts):
        # When the rule's params.sample_names is the wrong length the error
        # must say so explicitly. Without this check, downstream pandas
        # raises something cryptic about column-mismatch that doesn't point
        # at the rule definition mistake.
        counts_paths, _ = two_sample_counts
        scores_csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {"variants": ["p.(M1F)"], "mean": [-0.5], "sd": [0.1]}
        ).to_csv(scores_csv, index=False)
        with pytest.raises(ValueError) as excinfo:
            format_mavedb.format_scores(
                scores_path=str(scores_csv),
                backend="rosace",
                hgvs_col="variants",
                score_col="mean",
                sd_col="sd",
                output_path=str(tmp_path / "out_score.csv"),
                counts_paths=[str(p) for p in counts_paths],
                sample_names=["only_one"],
                counts_output_path=str(tmp_path / "out_count.csv"),
            )
        msg = str(excinfo.value)
        assert "sample_names" in msg
        # Both counts should appear in the error so the user can see
        # exactly which side has the wrong cardinality.
        assert "1" in msg and "2" in msg

    def test_scores_only_path_unchanged_when_no_counts(
        self, tmp_path, rosace_scores_csv
    ):
        # Belt-and-braces: counts_paths=None must take the original
        # scores-only emission path. No counts CSV emitted; score CSV has
        # exactly the input variants with no NaN-padding from union
        # expansion. Guarantees the new counts path doesn't accidentally
        # change behavior for callers that don't opt in.
        out = tmp_path / "out.csv"
        d = format_mavedb._BACKEND_DEFAULTS["rosace"]
        format_mavedb.format_scores(
            scores_path=str(rosace_scores_csv),
            backend="rosace",
            hgvs_col=d["hgvs_col"],
            score_col=d["score_col"],
            sd_col=d["sd_col"],
            output_path=str(out),
            counts_paths=None,
            sample_names=None,
            counts_output_path=None,
        )
        df = pd.read_csv(out)
        assert list(df.columns) == ["hgvs_pro", "score", "sd"]
        assert len(df) == 3
        assert not df["score"].isna().any()


# -----------------------------------------------------------------------------
# main() dispatch with counts wired. The rule now declares snakemake.input as
# a namespace with .scores + .counts, snakemake.output with .scores + .counts,
# and snakemake.params with backend + sample_names. main() must thread those
# through to format_scores so both CSVs land where the rule declared them.
# -----------------------------------------------------------------------------


class TestMainSnakemakeDispatchWithCounts:
    def test_main_emits_both_outputs_when_counts_wired(
        self, monkeypatch, tmp_path
    ):
        from types import SimpleNamespace

        scores_csv = tmp_path / "scores.csv"
        pd.DataFrame(
            {"variants": ["p.(M1F)"], "mean": [-0.5], "sd": [0.1]}
        ).to_csv(scores_csv, index=False)

        sample_csv = tmp_path / "A_R1_T0.csv"
        pd.DataFrame({"count": [42], "hgvs": ["p.(M1F)"]}).to_csv(
            sample_csv, index=False
        )

        score_out = tmp_path / "out_score.csv"
        count_out = tmp_path / "out_count.csv"

        mock_snakemake = SimpleNamespace(
            config={"mavedb": {}},
            params=SimpleNamespace(
                backend="rosace",
                sample_names=["A_R1_T0"],
            ),
            input=SimpleNamespace(
                scores=str(scores_csv),
                counts=[str(sample_csv)],
            ),
            output=SimpleNamespace(
                scores=str(score_out),
                counts=str(count_out),
            ),
        )
        monkeypatch.setattr(
            format_mavedb, "snakemake", mock_snakemake, raising=False
        )

        format_mavedb.main()

        score_df = pd.read_csv(score_out)
        count_df = pd.read_csv(count_out)
        # Single variant, present in both files, count column named after
        # the sample.
        assert score_df["hgvs_pro"].tolist() == ["p.Met1Phe"]
        assert count_df["hgvs_pro"].tolist() == ["p.Met1Phe"]
        assert count_df["A_R1_T0"].tolist() == [42]
