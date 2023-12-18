import pandas as pd
import logging
import snakemake

# The hierarchy of the enrich2 config file elements is as follows:
# experiment
#   conditions
#    tiles (if applicable: implemented as individual conditions)
#     replicates
#      timepoints/bins (individual samples)


def tiled_config(conditions, experiments, tsv_path, output_directory):
    """Generate the Enrich2 config file for tiled experiments."""
    enrich2_config = ["{", '\t"conditions": [']

    for condition in conditions:
        tiles = experiments.loc[(experiments["condition"] == condition)][
            "tile"
        ].unique()

        for tile in tiles:
            # open condition
            enrich2_config.extend(["\t{"])
            # condition name
            enrich2_config.extend(
                ['\t\t"name": "' + condition + "_tile" + str(tile) + '",']
            )
            # start selections (a.k.a. set of replicates)
            enrich2_config.extend(['\t\t"selections": ['])

            replicates = experiments.loc[
                (experiments["condition"] == condition) & (experiments["tile"] == tile)
            ]["replicate"].unique()

            for replicate in replicates:
                timepoints = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["replicate"] == replicate)
                    & (experiments["tile"] == tile)
                ]["time"].unique()

                # open library (a.k.a. single replicate)
                enrich2_config.extend(["\t\t{"])
                # start libraries
                enrich2_config.extend(['\t\t\t"libraries": ['])
                for time in timepoints:
                    sample_name = experiments.loc[
                        (experiments["condition"] == condition)
                        & (experiments["replicate"] == replicate)
                        & (experiments["time"] == time)
                        & (experiments["tile"] == tile)
                    ]["sample"].iloc[0]

                    # open file definition
                    name = "{}_rep{}_T{}_tile{}".format(
                        condition, str(replicate), str(time), str(tile)
                    )

                    enrich2_config.extend(["\t\t\t\t{"])

                    # close file definition
                    enrich2_config.extend(
                        [
                            '\t\t\t\t\t"counts file": "'
                            + tsv_path
                            + sample_name
                            + '.tsv",',
                            '\t\t\t\t\t"identifiers": {},',
                            '\t\t\t\t\t"name": "' + name + '",',
                            '\t\t\t\t\t"report filtered reads": false,',
                            '\t\t\t\t\t"timepoint": ' + str(time),
                        ]
                    )
                    if time == timepoints[-1]:
                        enrich2_config.extend(["\t\t\t\t}"])
                    else:
                        enrich2_config.extend(["\t\t\t\t},"])

                # closes and completes library definition
                enrich2_config.extend(["\t\t\t],"])
                enrich2_config.extend(
                    [
                        '\t\t\t"name": "'
                        + condition
                        + "_R"
                        + str(replicate)
                        + "_tile"
                        + str(tile)
                        + '"'
                    ]
                )
                if replicate == replicates[-1]:
                    enrich2_config.extend(["\t\t}"])
                else:
                    enrich2_config.extend(["\t\t},"])

                # closes replicate
            enrich2_config.extend(["\t\t]"])

            # closes condition. Checks for this being the last
            # condition in the file, as well.
            if condition == conditions[-1] and tile == tiles[-1]:
                enrich2_config.extend(["\t}"])
            else:
                enrich2_config.extend(["\t},"])

    enrich2_config.extend(["\t],"])
    enrich2_config.extend(['"name": "' + experiment_name + '",'])
    enrich2_config.extend(['"output directory": "' + output_directory + '"'])
    enrich2_config.extend(["}"])

    return enrich2_config


def untiled_config(conditions, experiments, tsv_path, output_directory):
    """Generate the Enrich2 config file for untiled experiments."""
    enrich2_config = ["{", '\t"conditions": [']

    for condition in conditions:
        # open condition
        enrich2_config.extend(["\t{"])
        # condition name
        enrich2_config.extend(['\t\t"name": "' + condition + '",'])
        # start selections (a.k.a. set of replicates)
        enrich2_config.extend(['\t\t"selections": ['])

        replicates = experiments.loc[(experiments["condition"] == condition)][
            "replicate"
        ].unique()

        for replicate in replicates:
            timepoints = experiments.loc[
                (experiments["condition"] == condition)
                & (experiments["replicate"] == replicate)
            ]["time"].unique()

            # open library (a.k.a. single replicate)
            enrich2_config.extend(["\t\t{"])
            # start libraries
            enrich2_config.extend(['\t\t\t"libraries": ['])
            for time in timepoints:
                sample_name = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["replicate"] == replicate)
                    & (experiments["time"] == time)
                ]["sample"].iloc[0]

                # open file definition
                name = "{}_rep{}_T{}".format(condition, str(replicate), str(time))

                enrich2_config.extend(["\t\t\t\t{"])

                # close file definition
                enrich2_config.extend(
                    [
                        '\t\t\t\t\t"counts file": "'
                        + tsv_path
                        + sample_name
                        + '.tsv",',
                        '\t\t\t\t\t"identifiers": {},',
                        '\t\t\t\t\t"name": "' + name + '",',
                        '\t\t\t\t\t"report filtered reads": false,',
                        '\t\t\t\t\t"timepoint": ' + str(time),
                    ]
                )
                if time == timepoints[-1]:
                    enrich2_config.extend(["\t\t\t\t}"])
                else:
                    enrich2_config.extend(["\t\t\t\t},"])

            # closes and completes library definition
            enrich2_config.extend(["\t\t\t],"])
            enrich2_config.extend(
                ['\t\t\t"name": "' + condition + "_R" + str(replicate) + '"']
            )
            if replicate == replicates[-1]:
                enrich2_config.extend(["\t\t}"])
            else:
                enrich2_config.extend(["\t\t},"])

            # closes replicate
        enrich2_config.extend(["\t\t]"])

        # closes condition
        if condition == conditions[-1]:
            enrich2_config.extend(["\t}"])
        else:
            enrich2_config.extend(["\t},"])

    enrich2_config.extend(["\t],"])
    enrich2_config.extend(['"name": "' + experiment_name + '",'])
    enrich2_config.extend(['"output directory": "' + output_directory + '"'])
    enrich2_config.extend(["}"])

    return enrich2_config


output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

experiment_name = snakemake.config["experiment"]
tsv_path = "results" + "/" + experiment_name + "/processed_counts/"

baseline_condition = snakemake.config["baseline_condition"]

output_directory = "results" + "/" + experiment_name + "/enrich/"

experiments = (
    pd.read_csv(snakemake.config["experiment_file"], header=0)
    .dropna(how="all")
    .set_index("sample", drop=False, verify_integrity=True)
)

conditions = experiments["condition"].unique().tolist()

# Do not generate scores for the baseline condition, if it exists.
if baseline_condition:
    try:
        conditions.remove(baseline_condition)
    except ValueError:
        logging.warning(
            "Baseline condition %s not found in experiment file.", baseline_condition
        )

# Generate the Enrich2 config file
if snakemake.config["tiled"]:
    enrich2_config_lines = tiled_config(conditions, experiments, tsv_path, output_directory)
else:
    enrich2_config_lines = untiled_config(conditions, experiments, tsv_path, output_directory)

with open(output_file, "w+") as f:
    for line in enrich2_config_lines:
        f.write(line + "\n")
