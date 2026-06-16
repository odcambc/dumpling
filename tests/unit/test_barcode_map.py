"""Unit tests for the barcode->variant map loader (barcode-counting mode).

Covers the decided behavior of script_utils.load_barcode_map: required-column
validation, barcode normalization, same-variant duplicate collapse, blank-row
dropping, BOM tolerance, and the empty-map guard. The conflict-resolution
policy (one tag -> multiple distinct variants) is exercised by the skipped test
at the bottom, which is unskipped once resolve_barcode_conflicts is implemented.
"""

import pandas as pd
import pytest

from workflow.rules.scripts.script_utils import (
    load_barcode_map,
    resolve_barcode_conflicts,
    validate_barcode_config,
)


def _write_map(tmp_path, rows, header="barcode,variant", name="bc_map.csv", encoding="utf-8"):
    """Write a tiny barcode-map CSV and return its path."""
    path = tmp_path / name
    lines = [header] + rows
    path.write_text("\n".join(lines) + "\n", encoding=encoding)
    return path


class TestLoadBarcodeMap:
    def test_simple_map_loads_to_dict(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)", "CCCC,p.(R11K)"])
        mapping, conflicts = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)", "CCCC": "p.(R11K)"}
        assert conflicts.empty

    def test_missing_barcode_column_raises(self, tmp_path):
        path = _write_map(tmp_path, ["p.(G10A)"], header="tag,variant")
        with pytest.raises(ValueError, match="missing required column"):
            load_barcode_map(path)

    def test_missing_variant_column_raises(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA"], header="barcode,hgvs")
        with pytest.raises(ValueError, match="missing required column"):
            load_barcode_map(path)

    def test_custom_variant_column_name(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)"], header="barcode,hgvs")
        mapping, _ = load_barcode_map(path, variant_column="hgvs")
        assert mapping == {"AAAA": "p.(G10A)"}

    def test_barcodes_are_uppercased_and_stripped(self, tmp_path):
        path = _write_map(tmp_path, ["  aAaA  ,p.(G10A)"])
        mapping, _ = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)"}

    def test_duplicate_same_variant_rows_collapse(self, tmp_path):
        # Same tag, same variant, listed twice -> one entry, no conflict.
        path = _write_map(tmp_path, ["AAAA,p.(G10A)", "AAAA,p.(G10A)"])
        mapping, conflicts = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)"}
        assert conflicts.empty

    def test_blank_variant_rows_are_dropped(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)", "CCCC,"])
        mapping, _ = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)"}

    def test_fully_blank_rows_ignored(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)", "", "  "])
        mapping, _ = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)"}

    def test_empty_map_raises(self, tmp_path):
        path = _write_map(tmp_path, ["CCCC,"])  # only an unusable row
        with pytest.raises(ValueError, match="no usable rows"):
            load_barcode_map(path)

    def test_bom_prefixed_map_loads(self, tmp_path):
        # Excel-exported CSVs carry a UTF-8 BOM; load_experiments tolerates it
        # via utf-8-sig and the barcode loader must too.
        path = _write_map(tmp_path, ["AAAA,p.(G10A)"], encoding="utf-8-sig")
        mapping, _ = load_barcode_map(path)
        assert mapping == {"AAAA": "p.(G10A)"}

    def test_ambiguous_barcode_dropped_and_reported(self, tmp_path):
        # AAAA names two distinct variants -> dropped (not in mapping), recorded
        # in the conflicts frame; the clean tag CCCC survives.
        path = _write_map(
            tmp_path, ["AAAA,p.(G10A)", "AAAA,p.(R11K)", "CCCC,p.(L12P)"]
        )
        mapping, conflicts = load_barcode_map(path)
        assert mapping == {"CCCC": "p.(L12P)"}
        assert list(conflicts["barcode"]) == ["AAAA"]
        assert conflicts.loc[0, "variants"] == "p.(G10A);p.(R11K)"

    def test_conflicts_written_to_duplicates_file(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)", "AAAA,p.(R11K)"])
        report = tmp_path / "nested" / "duped_barcodes.csv"
        mapping, conflicts = load_barcode_map(path, conflict_report_path=report)
        assert mapping == {}
        assert report.exists()  # parent dirs created as needed
        written = pd.read_csv(report)
        assert list(written["barcode"]) == ["AAAA"]

    def test_no_duplicates_file_when_no_conflicts(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA,p.(G10A)"])
        report = tmp_path / "duped_barcodes.csv"
        load_barcode_map(path, conflict_report_path=report)
        assert not report.exists()  # don't litter an empty report


class TestResolveBarcodeConflicts:
    def test_single_variant_is_kept(self):
        counts = pd.Series({"p.(G10A)": 3})
        assert resolve_barcode_conflicts("AAAA", counts) == "p.(G10A)"

    def test_multiple_variants_dropped(self):
        # drop-all policy: any tag naming >1 distinct variant is dropped.
        counts = pd.Series({"p.(G10A)": 9, "p.(R11K)": 1})
        assert resolve_barcode_conflicts("AAAA", counts) is None


class TestValidateBarcodeConfig:
    def _valid_map(self, tmp_path):
        return _write_map(tmp_path, ["AAAA,p.(G10A)"])

    def test_noop_when_not_barcoded(self, tmp_path):
        # No raise even with otherwise-broken barcode settings.
        validate_barcode_config({"barcoded": False, "barcode_pattern": "(", "barcode_start": 0})
        validate_barcode_config({})  # missing key entirely

    def test_missing_barcode_map_raises(self, tmp_path):
        with pytest.raises(ValueError, match="requires 'barcode_map'"):
            validate_barcode_config({"barcoded": True, "barcode_pattern": "^(.{4})"})

    def test_absent_barcode_map_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist"):
            validate_barcode_config(
                {
                    "barcoded": True,
                    "barcode_map": str(tmp_path / "nope.csv"),
                    "barcode_pattern": "^(.{4})",
                }
            )

    def test_map_missing_column_raises(self, tmp_path):
        path = _write_map(tmp_path, ["AAAA"], header="barcode,hgvs")
        with pytest.raises(ValueError, match="missing required column"):
            validate_barcode_config(
                {"barcoded": True, "barcode_map": str(path), "barcode_pattern": "^(.{4})"}
            )

    def test_both_extraction_specs_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="exactly one extraction spec"):
            validate_barcode_config(
                {
                    "barcoded": True,
                    "barcode_map": str(path),
                    "barcode_pattern": "^(.{4})",
                    "barcode_start": 0,
                    "barcode_length": 4,
                }
            )

    def test_neither_extraction_spec_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="exactly one extraction spec"):
            validate_barcode_config({"barcoded": True, "barcode_map": str(path)})

    def test_valid_pattern_passes(self, tmp_path):
        path = self._valid_map(tmp_path)
        validate_barcode_config(
            {"barcoded": True, "barcode_map": str(path), "barcode_pattern": "^(.{18})"}
        )

    def test_valid_fixed_position_passes(self, tmp_path):
        path = self._valid_map(tmp_path)
        validate_barcode_config(
            {
                "barcoded": True,
                "barcode_map": str(path),
                "barcode_start": 0,
                "barcode_length": 18,
            }
        )

    def test_pattern_without_capture_group_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="exactly one capture group"):
            validate_barcode_config(
                {"barcoded": True, "barcode_map": str(path), "barcode_pattern": "^.{18}"}
            )

    def test_pattern_with_two_capture_groups_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="exactly one capture group"):
            validate_barcode_config(
                {"barcoded": True, "barcode_map": str(path), "barcode_pattern": "(.{4})(.{4})"}
            )

    def test_invalid_regex_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="not a valid"):
            validate_barcode_config(
                {"barcoded": True, "barcode_map": str(path), "barcode_pattern": "^(.{4}"}
            )

    def test_negative_start_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="non-negative"):
            validate_barcode_config(
                {
                    "barcoded": True,
                    "barcode_map": str(path),
                    "barcode_start": -1,
                    "barcode_length": 4,
                }
            )

    def test_nonpositive_length_raises(self, tmp_path):
        path = self._valid_map(tmp_path)
        with pytest.raises(ValueError, match="must be positive"):
            validate_barcode_config(
                {
                    "barcoded": True,
                    "barcode_map": str(path),
                    "barcode_start": 0,
                    "barcode_length": 0,
                }
            )
