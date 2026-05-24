import json

import pandas as pd
import pytest

from workflow.rules.scripts.generate_enrich_configs import (
    remove_truncated_replicates,
    generate_config,
    remove_missing_t0,
    _run,
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
