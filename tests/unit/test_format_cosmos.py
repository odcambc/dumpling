# test_format_cosmos.py
#
# Tests for the cosmos (pimentellab/cosmos) score exporter. cosmos consumes a
# single wide CSV joining each phenotype's per-variant effect into
# beta_hat_N/se_hat_N columns; format_cosmos builds it from dumpling's
# per-condition score CSVs. See docs/cosmos_export_design.md.

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "workflow" / "rules" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import format_cosmos  # noqa: E402

# -----------------------------------------------------------------------------
# derive_cosmos_type: the one decision-bearing transform
# -----------------------------------------------------------------------------


class TestDeriveCosmosType:
    def test_substitutions_pass_through(self):
        # missense/synonymous/nonsense: the `type` label is already the cosmos
        # label; the substituted-AA `mutation` carries no length code.
        assert format_cosmos.derive_cosmos_type("missense", "F") == "missense"
        assert format_cosmos.derive_cosmos_type("synonymous", "S") == "synonymous"
        assert format_cosmos.derive_cosmos_type("nonsense", "X") == "nonsense"

    def test_synonymous_with_blank_mutation(self):
        # Synonymous rows may carry an empty/NaN mutation; must still pass through.
        assert format_cosmos.derive_cosmos_type("synonymous", "") == "synonymous"
        assert format_cosmos.derive_cosmos_type("synonymous", float("nan")) == (
            "synonymous"
        )
        assert format_cosmos.derive_cosmos_type("wt", None) == "wt"

    def test_insertion_fuses_length_no_separator(self):
        # "I_3" -> length 3, spelled "insertion3" (cosmos vignette convention).
        assert format_cosmos.derive_cosmos_type("insertion", "I_3") == "insertion3"
        assert format_cosmos.derive_cosmos_type("insertion", "I_1") == "insertion1"

    def test_deletion_fuses_length(self):
        assert format_cosmos.derive_cosmos_type("deletion", "D_2") == "deletion2"

    def test_insdel_keeps_distinct_label(self):
        # In-frame insdel (ID_<len>) is a distinct, length-resolved type, not
        # folded into insertion/deletion.
        assert format_cosmos.derive_cosmos_type("insdel", "ID_1") == "insdel1"

    def test_indel_detected_by_mutation_code_not_type_string(self):
        # Detection keys off the I_/D_/ID_ mutation prefix, so an oddly-labelled
        # `type` still gets length-resolved (uses whatever type label is given).
        assert format_cosmos.derive_cosmos_type("ins", "I_2") == "ins2"


# -----------------------------------------------------------------------------
# format_cosmos: multi-condition wide join
# -----------------------------------------------------------------------------


def _write_scores(path, rows):
    """rows: list of (variants, position, mutation, type, mean, sd)."""
    pd.DataFrame(
        rows, columns=["variants", "position", "mutation", "type", "mean", "sd"]
    ).to_csv(path, index=False)


ROSACE_COLS = format_cosmos._BACKEND_DEFAULTS["rosace"]


class TestFormatCosmos:
    def test_two_condition_outer_join(self, tmp_path):
        a = tmp_path / "cond_A_scores.csv"
        b = tmp_path / "cond_B_scores.csv"
        # V2/V3 shared; V1 only in A; V4 only in B.
        _write_scores(
            a,
            [
                ("p.(G2A)", 2, "A", "missense", 1.0, 0.1),
                ("p.(G3C)", 3, "C", "missense", 2.0, 0.2),
                ("p.(M1del)", 1, "D_1", "deletion", 9.0, 0.9),
            ],
        )
        _write_scores(
            b,
            [
                ("p.(G2A)", 2, "A", "missense", -1.0, 0.3),
                ("p.(G3C)", 3, "C", "missense", -2.0, 0.4),
                ("p.(K5R)", 5, "R", "missense", 5.0, 0.5),
            ],
        )
        out = tmp_path / "exp_cosmos.csv"
        format_cosmos.format_cosmos([str(a), str(b)], "rosace", ROSACE_COLS, str(out))

        df = pd.read_csv(out)
        # Required cosmos columns, in order.
        assert list(df.columns) == [
            "variants",
            "group",
            "type",
            "beta_hat_1",
            "se_hat_1",
            "beta_hat_2",
            "se_hat_2",
        ]
        # Outer join: union of all variants.
        assert set(df["variants"]) == {"p.(G2A)", "p.(G3C)", "p.(M1del)", "p.(K5R)"}

        g2a = df[df["variants"] == "p.(G2A)"].iloc[0]
        assert g2a["beta_hat_1"] == 1.0 and g2a["beta_hat_2"] == -1.0
        assert g2a["group"] == 2 and g2a["type"] == "missense"

        # V1 only in cond A (slot 1): beta_hat_2 is blank; deletion length fused.
        m1 = df[df["variants"] == "p.(M1del)"].iloc[0]
        assert m1["beta_hat_1"] == 9.0 and pd.isna(m1["beta_hat_2"])
        assert m1["type"] == "deletion1"

        # V4 only in cond B (slot 2): beta_hat_1 blank, group/type coalesced from B.
        k5 = df[df["variants"] == "p.(K5R)"].iloc[0]
        assert pd.isna(k5["beta_hat_1"]) and k5["beta_hat_2"] == 5.0
        assert k5["group"] == 5

    def test_slot_order_follows_input_order(self, tmp_path):
        a = tmp_path / "a.csv"
        b = tmp_path / "b.csv"
        _write_scores(a, [("p.(G2A)", 2, "A", "missense", 1.0, 0.1)])
        _write_scores(b, [("p.(G2A)", 2, "A", "missense", 2.0, 0.2)])
        out = tmp_path / "o.csv"
        # b first -> b is beta_hat_1.
        format_cosmos.format_cosmos([str(b), str(a)], "rosace", ROSACE_COLS, str(out))
        row = pd.read_csv(out).iloc[0]
        assert row["beta_hat_1"] == 2.0 and row["beta_hat_2"] == 1.0

    def test_duplicate_variant_in_condition_raises(self, tmp_path):
        a = tmp_path / "a.csv"
        _write_scores(
            a,
            [
                ("p.(G2A)", 2, "A", "missense", 1.0, 0.1),
                ("p.(G2A)", 2, "A", "missense", 1.5, 0.1),
            ],
        )
        with pytest.raises(ValueError, match="duplicate variants"):
            format_cosmos.format_cosmos(
                [str(a)], "rosace", ROSACE_COLS, str(tmp_path / "o.csv")
            )

    def test_missing_column_raises_named_error(self, tmp_path):
        a = tmp_path / "a.csv"
        pd.DataFrame(
            [("p.(G2A)", 2, "missense", 1.0)],
            columns=["variants", "position", "type", "mean"],  # no sd, no mutation
        ).to_csv(a, index=False)
        with pytest.raises(ValueError, match="Expected columns not found"):
            format_cosmos.format_cosmos(
                [str(a)], "rosace", ROSACE_COLS, str(tmp_path / "o.csv")
            )

    def test_lilace_backend_column_names(self, tmp_path):
        cols = format_cosmos._BACKEND_DEFAULTS["lilace"]
        a = tmp_path / "lilace.csv"
        pd.DataFrame(
            [("p.(G2A)", 2, "A", "missense", 1.0, 0.1)],
            columns=["variant", "position", "mutation", "type", "effect", "effect_se"],
        ).to_csv(a, index=False)
        out = tmp_path / "o.csv"
        format_cosmos.format_cosmos([str(a)], "lilace", cols, str(out))
        df = pd.read_csv(out)
        assert df.iloc[0]["beta_hat_1"] == 1.0 and df.iloc[0]["se_hat_1"] == 0.1
