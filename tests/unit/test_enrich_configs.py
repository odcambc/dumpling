import pandas as pd
from workflow.rules.scripts.generate_enrich_configs import (
    remove_truncated_replicates,
    generate_config,
    remove_missing_t0,
)


def test_remove_truncated_replicates():
    data = {
        "condition": ["A", "A", "A", "B", "B"],
        "tile": [1, 1, 1, 2, 2],
        "replicate": [1, 1, 1, 1, 1],
        "time": [0, 1, 2, 0, 1],  # Condition B has fewer than 2 timepoints
        "sample": ["s1", "s2", "s3", "s4", "s5"],
    }
    df = pd.DataFrame(data)
    conditions = ["A", "B"]

    filtered_df = remove_truncated_replicates(df, conditions, tiled=False)

    assert len(filtered_df) == 3  # Only condition A should remain
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
