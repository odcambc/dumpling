import pandas as pd
import logging

import yaml


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

experiments = pd.read_csv(snakemake.config["experiment_file"], header=0).set_index(
    "sample", drop=False, verify_integrity=True
)

enrich2_config = ["{", '\t"conditions": [']

conditions = experiments["condition"].unique()
replicates = experiments["replicate"].unique()

# Do not generate scores for the baseline condition, if it exists.
if baseline_condition:
    try:
        conditions = conditions.remove(baseline_condition)
    except ValueError:
        logging.warning("Baseline condition %s not found in experiment file.", baseline_condition)

for condition in conditions:
    # open condition

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
            ]["sample"][0]

            # open file definition
            name = "{}_rep{}_T{}".format(condition, str(replicate), str(time))

            enrich2_config.extend(["\t\t\t\t{"])

            # close file definition
            enrich2_config.extend(
                [
                    '\t\t\t\t\t"counts file": "' + tsv_path + sample_name + '.tsv",',
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

with open(output_file, "w+") as f:
    for line in enrich2_config:
        f.write(line + "\n")