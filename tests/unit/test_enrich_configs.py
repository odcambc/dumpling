import json

import pandas as pd
import pytest

from workflow.rules.scripts.generate_enrich_configs import (
    _run,
    enrich_library_label,
    enrich_selection_label,
    expected_enrich_h5_basenames,
    generate_config,
    remove_missing_t0,
    remove_truncated_replicates,
)


def test_remove_truncated_replicates_keeps_two_timepoint_replicates():
    """Enrich2 needs T0 + at least one later sample (2 timepoints total) to
    compute a ratio. A replicate with exactly 2 timepoints is valid and must
    be preserved — the prior `<= 2` check removed it incorrectly."""
    data = {
        "condition": ["A", "A", "A", "B", "B"],
        "tile": [1, 1, 1, 2, 2],
        "replicate": [1, 1, 1, 1, 1],
        "time": [0, 1, 2, 0, 1],  # B has 2 timepoints — valid.
        "sample": ["s1", "s2", "s3", "s4", "s5"],
    }
    df = pd.DataFrame(data)
    conditions = ["A", "B"]

    filtered_df = remove_truncated_replicates(df, conditions, tiled=False)

    # Both A and B survive: A has 3 timepoints, B has 2.
    assert len(filtered_df) == 5
    assert set(filtered_df["condition"].unique()) == {"A", "B"}


def test_remove_truncated_replicates_drops_single_timepoint_replicates():
    """A replicate with only 1 timepoint can't form a ratio and must be
    dropped."""
    data = {
        "condition": ["A", "A", "A", "B"],
        "tile": [1, 1, 1, 2],
        "replicate": [1, 1, 1, 1],
        "time": [0, 1, 2, 0],  # B has only T0 — truly truncated.
        "sample": ["s1", "s2", "s3", "s4"],
    }
    df = pd.DataFrame(data)
    conditions = ["A", "B"]

    filtered_df = remove_truncated_replicates(df, conditions, tiled=False)

    assert len(filtered_df) == 3
    assert "B" not in filtered_df["condition"].unique()


def test_remove_missing_t0():
    data = {
        "condition": ["A", "A", "B", "B"],
        "tile": [1, 1, 2, 2],
        "replicate": [1, 1, 1, 1],
        "time": [1, 2, 0, 1],  # Condition A has no T0
        "sample": ["s1", "s2", "s3", "s4"],
    }
    df = pd.DataFrame(data)
    conditions = ["A", "B"]

    filtered_df = remove_missing_t0(df, conditions, tiled=False)

    assert len(filtered_df) == 2  # Only condition B should remain
    assert "A" not in filtered_df["condition"].unique()


def test_generate_config_tiled():
    data = {
        "condition": ["A", "A"],
        "tile": [1, 1],
        "replicate": [1, 1],
        "time": [0, 1],
        "sample": ["s1", "s2"],
    }
    df = pd.DataFrame(data)
    conditions = ["A"]
    tsv_path = "/path/to/tsv/"
    output_directory = "/output/directory/"
    experiment_name = "test_experiment"

    config = generate_config(
        conditions,
        df,
        tsv_path,
        output_directory,
        tiled=True,
        experiment_name=experiment_name,
    )

    assert config[0] == "{"
    assert any('"name": "A_tile1"' in line for line in config)
    assert any(f'"counts file": "{tsv_path}s1.tsv"' in line for line in config)
    assert any(f'"output directory": "{output_directory}"' in line for line in config)


def test_generate_config_untiled():
    data = {
        "condition": ["A", "A"],
        "tile": [None, None],
        "replicate": [1, 1],
        "time": [0, 1],
        "sample": ["s1", "s2"],
    }
    df = pd.DataFrame(data)
    conditions = ["A"]
    tsv_path = "/path/to/tsv/"
    output_directory = "/output/directory/"
    experiment_name = "test_experiment"

    config = generate_config(
        conditions,
        df,
        tsv_path,
        output_directory,
        tiled=False,
        experiment_name=experiment_name,
    )

    assert config[0] == "{"
    assert any('"name": "A"' in line for line in config)
    assert any(f'"counts file": "{tsv_path}s1.tsv"' in line for line in config)
    assert any(f'"output directory": "{output_directory}"' in line for line in config)


# ---------------------------------------------------------------------------
# Empty-condition handling end-to-end via _run.
#
# remove_truncated_replicates / remove_missing_t0 can drop every replicate of
# a condition. Before this fix, _run still passed the *original* conditions
# list into generate_config, so the dropped condition would appear as a
# stanza with an empty "selections": [] — and the trailing-comma logic could
# emit invalid JSON if the empty condition happened to be last.
# ---------------------------------------------------------------------------


def _write_experiment_csv(tmp_path, rows):
    csv = tmp_path / "experiment.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    return csv


def _make_snakemake(mock_snakemake, experiment_csv, output_path, *, remove_zeros=False):
    return mock_snakemake(
        config={
            "experiment": "test_exp",
            "experiment_file": str(experiment_csv),
            "tiled": False,
            "baseline_condition": "",
        },
        params={"remove_zeros": remove_zeros},
        output=[str(output_path)],
        log=["/dev/null"],
    )


def test_run_drops_condition_filtered_out_by_replicate_pruning(tmp_path, mock_snakemake):
    """Condition B's only replicate has <2 timepoints; after filtering it
    must not appear in the generated config, and the surviving JSON must
    parse cleanly."""
    rows = [
        # Condition A: valid (T0, T1, T2)
        {"sample": "A_T0", "condition": "A", "replicate": 1, "time": 0},
        {"sample": "A_T1", "condition": "A", "replicate": 1, "time": 1},
        {"sample": "A_T2", "condition": "A", "replicate": 1, "time": 2},
        # Condition B: only one timepoint — gets dropped by remove_truncated_replicates
        {"sample": "B_T0", "condition": "B", "replicate": 1, "time": 0},
    ]
    csv = _write_experiment_csv(tmp_path, rows)
    out = tmp_path / "config.json"

    _run(_make_snakemake(mock_snakemake, csv, out))

    text = out.read_text()
    parsed = json.loads(text)  # must be valid JSON

    names = {c["name"] for c in parsed["conditions"]}
    assert names == {"A"}, f"Condition B should have been dropped; got {names}"
    # And the surviving condition has a non-empty selections list.
    assert len(parsed["conditions"][0]["selections"]) > 0


def test_run_raises_when_all_conditions_dropped(tmp_path, mock_snakemake):
    """If filtering removes every replicate of every condition, _run must
    fail loudly rather than emit an empty `conditions: []` config that
    Enrich2 will choke on later."""
    rows = [
        # Both conditions have only a single timepoint — both get dropped.
        {"sample": "A_T0", "condition": "A", "replicate": 1, "time": 0},
        {"sample": "B_T0", "condition": "B", "replicate": 1, "time": 0},
    ]
    csv = _write_experiment_csv(tmp_path, rows)
    out = tmp_path / "config.json"

    with pytest.raises(ValueError, match="No conditions remain"):
        _run(_make_snakemake(mock_snakemake, csv, out))


# ---------------------------------------------------------------------------
# expected_enrich_h5_basenames: predicts the .h5 stores Enrich2 writes so
# run_enrich can declare them as temp() outputs (issue #16). Must match what
# generate_config emits as object names, hence the shared label helpers.
# ---------------------------------------------------------------------------


class TestEnrichLabels:
    def test_selection_label_untiled(self):
        assert enrich_selection_label("cond_A", 1, use_tile=False) == "cond_A_R1"

    def test_selection_label_tiled(self):
        assert (
            enrich_selection_label("cond_A", 1, use_tile=True, tile=2)
            == "cond_A_R1_tile2"
        )

    def test_library_label_untiled(self):
        assert (
            enrich_library_label("cond_A", 1, 0, use_tile=False) == "cond_A_rep1_T0"
        )

    def test_library_label_tiled(self):
        assert (
            enrich_library_label("cond_A", 1, 0, use_tile=True, tile=2)
            == "cond_A_rep1_T0_tile2"
        )


class TestExpectedEnrichH5Basenames:
    def test_matches_committed_example_layout(self):
        """Pinned against the real results/example_experiment/enrich tree:
        cond_A has rep1 (T0-T3) and rep2 (T0-T2); cond_B has rep1 (T0-T3).
        The predicted store set must equal exactly the .h5 files Enrich2 wrote
        there — this is what guards the temp() declaration against drift."""
        rows = []
        for time in (0, 1, 2, 3):
            rows.append({"condition": "cond_A", "replicate": 1, "time": time})
        for time in (0, 1, 2):
            rows.append({"condition": "cond_A", "replicate": 2, "time": time})
        for time in (0, 1, 2, 3):
            rows.append({"condition": "cond_B", "replicate": 1, "time": time})
        for i, row in enumerate(rows):
            row["sample"] = f"s{i}"
        df = pd.DataFrame(rows)

        names = expected_enrich_h5_basenames(
            df, ["cond_A", "cond_B"], tiled=False, experiment_name="example_experiment"
        )

        assert set(names) == {
            "example_experiment_exp.h5",
            "cond_A_R1_sel.h5",
            "cond_A_R2_sel.h5",
            "cond_A_rep1_T0_lib.h5",
            "cond_A_rep1_T1_lib.h5",
            "cond_A_rep1_T2_lib.h5",
            "cond_A_rep1_T3_lib.h5",
            "cond_A_rep2_T0_lib.h5",
            "cond_A_rep2_T1_lib.h5",
            "cond_A_rep2_T2_lib.h5",
            "cond_B_R1_sel.h5",
            "cond_B_rep1_T0_lib.h5",
            "cond_B_rep1_T1_lib.h5",
            "cond_B_rep1_T2_lib.h5",
            "cond_B_rep1_T3_lib.h5",
        }

    def test_excludes_filtered_out_replicates(self):
        """A replicate dropped by the T0 / minimum-timepoint filtering must not
        contribute store files — otherwise temp() would declare an output
        Enrich2 never writes and the run would fail."""
        rows = [
            # cond_A rep1: valid (T0, T1)
            {"sample": "a0", "condition": "cond_A", "replicate": 1, "time": 0},
            {"sample": "a1", "condition": "cond_A", "replicate": 1, "time": 1},
            # cond_A rep2: only one timepoint -> dropped by truncation filter
            {"sample": "a2", "condition": "cond_A", "replicate": 2, "time": 0},
            # cond_B rep1: no T0 -> dropped by missing-T0 filter
            {"sample": "b1", "condition": "cond_B", "replicate": 1, "time": 1},
            {"sample": "b2", "condition": "cond_B", "replicate": 1, "time": 2},
        ]
        df = pd.DataFrame(rows)
        names = expected_enrich_h5_basenames(
            df, ["cond_A", "cond_B"], tiled=False, experiment_name="exp"
        )
        assert set(names) == {
            "exp_exp.h5",
            "cond_A_R1_sel.h5",
            "cond_A_rep1_T0_lib.h5",
            "cond_A_rep1_T1_lib.h5",
        }

    def test_agrees_with_generated_config_object_names(self):
        """Anti-drift: every selection/library object name generate_config puts
        in the config JSON must have a matching _sel.h5 / _lib.h5 in the
        predicted store set, since Enrich2 names each store after its object."""
        rows = [
            {"sample": "a0", "condition": "cond_A", "replicate": 1, "time": 0},
            {"sample": "a1", "condition": "cond_A", "replicate": 1, "time": 1},
            {"sample": "a2", "condition": "cond_A", "replicate": 1, "time": 2},
        ]
        df = pd.DataFrame(rows)
        config_lines = generate_config(
            ["cond_A"], df, "/tsv/", "/out/", tiled=False, experiment_name="exp"
        )
        config_names = {
            line.split('"name": "')[1].rstrip('",').rstrip('"')
            for line in config_lines
            if '"name": "' in line
        }
        stores = set(
            expected_enrich_h5_basenames(
                df, ["cond_A"], tiled=False, experiment_name="exp"
            )
        )
        # Selection object "cond_A_R1" -> "cond_A_R1_sel.h5";
        # library object "cond_A_rep1_T0" -> "cond_A_rep1_T0_lib.h5".
        for name in config_names:
            if name.startswith("cond_A_R") and "_rep" not in name:
                assert f"{name}_sel.h5" in stores, name
            elif "_rep" in name:
                assert f"{name}_lib.h5" in stores, name
