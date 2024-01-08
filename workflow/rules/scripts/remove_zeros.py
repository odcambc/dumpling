import logging
import pathlib

import pandas as pd

# Set up logging
log_file = snakemake.log[0]

if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

logging.debug("Performing unobserved variant removal.")


experiment_name = snakemake.config["experiment"]
output_dir = "results/" + experiment_name + "/processed_counts/removed_zeros/"
input_dir = "results/" + experiment_name + "/processed_counts/"


def remove_zeros(experiment_list, output_dir):
    """Remove variants that are not observed in any of the experiments."""

    enrich_file_list = [input_dir + s + ".tsv" for s in experiment_list]

    unobserved = remove_zeros_enrich(enrich_file_list, output_dir)

    unobserved_file = (
        output_dir + "rejected/" + experiment_name + "_unobserved_variants.csv"
    )

    p = pathlib.Path(unobserved_file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        logging.debug("Writing unobserved variants to %s.", p)
        for variant in unobserved:
            f.write("%s\n" % variant)



def remove_zeros_enrich(enrich_file_list, output_dir):
    # Input is a list of files containing Enrich2 formatted tsvs.
    # This function combines them all into a single df, then removes the
    # variants that are missing across all experiments, and rewrites the tsvs.

    df_list = [
        pd.read_csv(f, sep="\t", names=["hgvs", f], header=0).set_index("hgvs")
        for f in enrich_file_list
    ]

    combined_enrich_df = pd.concat(df_list, join="outer", axis=1)

    unobserved_variants = combined_enrich_df[
        (combined_enrich_df == 0).all(axis=1)
    ].index.to_list()

    combined_enrich_df = combined_enrich_df[(combined_enrich_df != 0).any(axis=1)]

    # Write out each file to the zero-filtered directory

    for enrich_file in enrich_file_list:
        file_name = pathlib.Path(enrich_file).name
        p = pathlib.Path(output_dir, file_name)
        p.parent.mkdir(parents=True, exist_ok=True)

        with p.open("w+") as f:
            combined_enrich_df.to_csv(
                f, columns=[enrich_file], header=["count"], index=True, sep="\t"
            )
            logging.debug("Writing zero-filtered file %s to %s.", enrich_file, p)


    return unobserved_variants


# Process the experiments according to their relationships
with open(snakemake.config["experiment_file"]) as f:
    experiments = (
        pd.read_csv(f, header=0)
        .dropna(how="all")
        .set_index("sample", drop=False, verify_integrity=True)
    )

for condition in experiments["condition"].unique():
    logging.debug("Condition: %s", condition)
    tile_list = experiments.loc[(experiments["condition"] == condition)]["tile"]
    for tile in tile_list.unique():
        logging.debug("Tile: %s", tile)
        replicate_list = experiments.loc[
            (experiments["condition"] == condition) & (experiments["tile"] == tile)
        ]["replicate"]
        for replicate in replicate_list.unique():
            logging.debug("Replicate: %s", replicate)
            experiment_list = experiments.loc[
                (experiments["condition"] == condition)
                & (experiments["replicate"] == replicate)
                & (experiments["tile"] == tile)
            ]["sample"].tolist()
            file_list = experiments.loc[
                (experiments["condition"] == condition)
                & (experiments["replicate"] == replicate)
                & (experiments["tile"] == tile)
            ]["file"].tolist()
            logging.debug("Samples: ")
            logging.debug(experiment_list)

            remove_zeros(experiment_list, output_dir)
