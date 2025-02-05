import pandas as pd
import logging
from snakemake.script import snakemake

# The hierarchy of the enrich2 config file elements is as follows:
# experiment
#   conditions
#    tiles (if applicable: implemented as individual conditions)
#     replicates
#      timepoints/bins (individual samples)


def remove_truncated_replicates(experiments, conditions, tiled):
    """
    Remove replicates with fewer than two timepoints/bins.

    Accepts as input a dataframe of experiment metadata.
    Returns a dataframe of experiment metadata with all replicates
    that have two or fewer timepoints removed.
    """
    has_tile_column = "tile" in experiments.columns

    for condition in conditions:
        # Determine the tile values we loop over
        if tiled and has_tile_column:
            tiles = experiments.loc[
                experiments["condition"] == condition, "tile"
            ].unique()
        else:
            # Treat as if there's only one tile when 'tile' doesn't exist or it's not tiled
            tiles = [None]

        for tile in tiles:
            # Subset for a given condition/tile (or no tile)
            if tile is None or not (tiled and has_tile_column):
                condition_subset = experiments.loc[
                    experiments["condition"] == condition
                ]
            else:
                condition_subset = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["tile"] == tile)
                ]

            replicates = condition_subset["replicate"].unique()

            for replicate in replicates:
                rep_subset = condition_subset.loc[
                    condition_subset["replicate"] == replicate
                ]
                timepoints = rep_subset["time"].unique()

                if len(timepoints) <= 2:
                    # Remove all rows matching this condition+replicate (+tile if present)
                    if tile is None or not (tiled and has_tile_column):
                        experiments = experiments.loc[
                            (experiments["condition"] != condition)
                            | (experiments["replicate"] != replicate)
                        ]
                    else:
                        experiments = experiments.loc[
                            (experiments["condition"] != condition)
                            | (experiments["replicate"] != replicate)
                            | (experiments["tile"] != tile)
                        ]
                    logging.warning(
                        "Replicate %s in condition %s has fewer than two timepoints. "
                        "Removing from analysis.",
                        replicate,
                        condition,
                    )

    return experiments


def remove_missing_t0(experiments, conditions, tiled):
    """
    Remove replicates without a T0 sample.

    Accepts as input:
      - experiments: DataFrame of experiment metadata.
      - conditions: List of condition names.
      - tiled: Boolean indicating whether the experiment is tiled or not.

    Returns:
      experiments: DataFrame with all replicates lacking T0 removed.
    """
    has_tile_column = "tile" in experiments.columns

    for condition in conditions:
        if tiled and has_tile_column:
            tiles = experiments.loc[
                experiments["condition"] == condition, "tile"
            ].unique()
        else:
            tiles = [None]

        for tile in tiles:
            if tile is None or not (tiled and has_tile_column):
                condition_subset = experiments.loc[
                    experiments["condition"] == condition
                ]
            else:
                condition_subset = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["tile"] == tile)
                ]

            replicates = condition_subset["replicate"].unique()

            for replicate in replicates:
                rep_subset = condition_subset.loc[
                    condition_subset["replicate"] == replicate
                ]
                timepoints = rep_subset["time"].unique()

                if 0 not in timepoints:
                    if tile is None or not (tiled and has_tile_column):
                        experiments = experiments.loc[
                            (experiments["condition"] != condition)
                            or (experiments["replicate"] != replicate)
                        ]
                    else:
                        experiments = experiments.loc[
                            (experiments["condition"] != condition)
                            | (experiments["replicate"] != replicate)
                            | (experiments["tile"] != tile)
                        ]
                    logging.warning(
                        "Replicate %s in condition %s has no T0 sample. "
                        "Removing from analysis.",
                        replicate,
                        condition,
                    )

    return experiments


def generate_config(
    conditions, experiments, tsv_path, output_directory, experiment_name, tiled
):
    """Generate the Enrich2 config file for experiments."""
    enrich2_config = ["{", '\t"conditions": [']

    for condition in conditions:
        tiles = (
            experiments.loc[experiments["condition"] == condition, "tile"].unique()
            if (tiled and "tile" in experiments.columns)
            else [None]
        )

        for tile in tiles:
            enrich2_config.append("\t{")
            condition_name = (
                f'"name": "{condition}_tile{tile}"'
                if (tiled and "tile" in experiments.columns)
                else f'"name": "{condition}"'
            )
            enrich2_config.append(f"\t\t{condition_name},")
            enrich2_config.append('\t\t"selections": [')

            replicates = experiments.loc[
                (experiments["condition"] == condition)
                & (
                    (experiments["tile"] == tile)
                    if (tiled and "tile" in experiments.columns)
                    else True
                ),
                "replicate",
            ].unique()

            for replicate in replicates:
                timepoints = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["replicate"] == replicate)
                    & (
                        (experiments["tile"] == tile)
                        if (tiled and "tile" in experiments.columns)
                        else True
                    ),
                    "time",
                ].unique()

                enrich2_config.append("\t\t{")
                enrich2_config.append('\t\t\t"libraries": [')

                for time in timepoints:
                    sample_name = experiments.loc[
                        (experiments["condition"] == condition)
                        & (experiments["replicate"] == replicate)
                        & (experiments["time"] == time)
                        & (
                            (experiments["tile"] == tile)
                            if (tiled and "tile" in experiments.columns)
                            else True
                        ),
                        "sample",
                    ].iloc[0]
                    name = (
                        f"{condition}_rep{replicate}_T{time}_tile{tile}"
                        if (tiled and "tile" in experiments.columns)
                        else f"{condition}_rep{replicate}_T{time}"
                    )

                    enrich2_config.append("\t\t\t\t{")
                    enrich2_config.append(
                        f'\t\t\t\t\t"counts file": "{tsv_path}{sample_name}.tsv",'
                    )
                    enrich2_config.append('\t\t\t\t\t"identifiers": {},')
                    enrich2_config.append(f'\t\t\t\t\t"name": "{name}",')
                    enrich2_config.append('\t\t\t\t\t"report filtered reads": false,')
                    enrich2_config.append(f'\t\t\t\t\t"timepoint": {time}')
                    enrich2_config.append(
                        "\t\t\t\t}," if time != timepoints[-1] else "\t\t\t\t}"
                    )

                enrich2_config.append("\t\t\t],")
                enrich2_config.append(
                    f'\t\t\t"name": "{condition}_R{replicate}_tile{tile}"'
                    if (tiled and "tile" in experiments.columns)
                    else f'\t\t\t"name": "{condition}_R{replicate}"'
                )
                enrich2_config.append(
                    "\t\t}," if replicate != replicates[-1] else "\t\t}"
                )

            enrich2_config.append("\t\t]")
            enrich2_config.append(
                "\t},"
                if (
                    condition != conditions[-1]
                    or (tiled and "tile" in experiments.columns and tile != tiles[-1])
                )
                else "\t}"
            )

    enrich2_config.append("\t],")
    enrich2_config.append(f'"name": "{experiment_name}",')
    enrich2_config.append(f'"output directory": "{output_directory}"')
    enrich2_config.append("}")

    return enrich2_config


def main():
    """
    Handle the main logic of the script.
    """
    # Set up logging
    log_file = snakemake.log[0]
    if log_file:
        logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG)

    experiment_name = snakemake.config["experiment"]
    tiled = snakemake.config["tiled"]

    # Decide which directory path to use
    if snakemake.params["remove_zeros"]:
        tsv_path = f"results/{experiment_name}/processed_counts/removed_zeros/"
    else:
        tsv_path = f"results/{experiment_name}/processed_counts/"

    baseline_condition = snakemake.config["baseline_condition"]
    output_directory = f"results/{experiment_name}/enrich/"
    output_file = snakemake.output[0]

    # Read in the experiments file
    experiments = (
        pd.read_csv(snakemake.config["experiment_file"], header=0)
        .dropna(how="all")
        .set_index("sample", drop=False, verify_integrity=True)
    )

    # List of condition names, ignoring baseline condition if present
    global conditions
    conditions = experiments["condition"].unique().tolist()
    if baseline_condition:
        try:
            conditions.remove(baseline_condition)
        except ValueError:
            logging.warning(
                "Baseline condition %s not found in experiment file.",
                baseline_condition,
            )

    # First remove replicates with no T0 sample, then replicates with <2 timepoints
    experiments_filtered = remove_missing_t0(experiments, conditions, tiled)
    experiments_filtered = remove_truncated_replicates(
        experiments_filtered, conditions, tiled
    )

    # Generate the Enrich2 config
    enrich2_config_lines = generate_config(
        conditions,
        experiments_filtered,
        tsv_path,
        output_directory,
        experiment_name,
        tiled,
    )

    # Write out the generated config file
    with open(output_file, "w+") as f:
        f.write("\n".join(enrich2_config_lines))

    logging.info("Enrich2 config file written to %s.", output_file)


if __name__ == "__main__":
    main()
