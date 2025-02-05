import logging
from pathlib import Path

import pandas as pd
from snakemake.script import snakemake


def remove_zeros_enrich(enrich_file_list, output_dir):
    """
    Combine Enrich2 TSV files, remove unobserved variants across all experiments,
    and rewrite the filtered TSV files.

    :param enrich_file_list: List of input file paths in Enrich2 format
    :param output_dir: Directory to place filtered output
    :return: List of unobserved variant IDs
    """
    # Build a combined DataFrame from all inputs
    df_list = [
        pd.read_csv(f, sep="\t", names=["hgvs", f], header=0).set_index("hgvs")
        for f in enrich_file_list
    ]
    combined_enrich_df = pd.concat(df_list, join="outer", axis=1)

    # Identify unobserved variants (those with all zeros)
    unobserved_variants = combined_enrich_df[
        (combined_enrich_df == 0).all(axis=1)
    ].index.to_list()

    # Keep only variants that are observed in at least one file
    combined_enrich_df = combined_enrich_df[(combined_enrich_df != 0).any(axis=1)]

    # Write out each file with filtered data
    for enrich_file in enrich_file_list:
        file_name = Path(enrich_file).name
        output_path = Path(output_dir, file_name)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write the single-column file for that sample
        combined_enrich_df[[enrich_file]].to_csv(
            output_path, header=["count"], index=True, sep="\t"
        )
        logging.debug("Writing zero-filtered file %s to %s.", enrich_file, output_path)

    return unobserved_variants


def remove_zeros(experiment_list, input_dir, output_dir, experiment_name):
    """
    Remove variants that are not observed in any of the experiments in `experiment_list`.
    Writes the list of unobserved variants to a separate file.
    """
    enrich_file_list = [f"{input_dir}{s}.tsv" for s in experiment_list]
    unobserved = remove_zeros_enrich(enrich_file_list, output_dir)

    unobserved_file = Path(
        output_dir, "rejected", f"{experiment_name}_unobserved_variants.csv"
    )
    unobserved_file.parent.mkdir(parents=True, exist_ok=True)

    with unobserved_file.open("w+") as f:
        logging.debug("Writing unobserved variants to %s.", unobserved_file)
        f.writelines(f"{variant}\n" for variant in unobserved)


def main():
    """
    Main function to perform unobserved variant removal for each condition/tile/replicate grouping.
    - Sets up logging
    - Reads config from Snakemake
    - Checks if 'tile' column exists for optional tiling
    - Calls remove_zeros() for each replicate group
    """
    # Set up logging
    log_file = snakemake.log[0]
    if log_file:
        logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG)

    logging.debug("Performing unobserved variant removal.")

    experiment_name = snakemake.config["experiment"]
    experiment_file = snakemake.config["experiment_file"]
    input_dir = f"results/{experiment_name}/processed_counts/enrich_format/"
    output_dir = f"results/{experiment_name}/processed_counts/removed_zeros/"

    # Read experiments
    experiments = (
        pd.read_csv(experiment_file, header=0)
        .dropna(how="all")
        .set_index("sample", drop=False, verify_integrity=True)
    )

    # Check if 'tile' column exists
    has_tile = "tile" in experiments.columns

    # Process experiments according to condition/tile/replicate
    for condition in experiments["condition"].unique():
        logging.debug("Condition: %s", condition)

        if has_tile:
            tile_values = experiments.loc[
                experiments["condition"] == condition, "tile"
            ].unique()
        else:
            # If no tile column, treat as one group
            tile_values = [None]

        for tile in tile_values:
            logging.debug("Tile: %s", tile if tile is not None else "no tile")
            if tile is None or not has_tile:
                condition_subset = experiments.loc[
                    experiments["condition"] == condition
                ]
            else:
                condition_subset = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["tile"] == tile)
                ]

            for replicate in condition_subset["replicate"].unique():
                logging.debug("Replicate: %s", replicate)
                # Filter by replicate
                replicate_subset = condition_subset.loc[
                    condition_subset["replicate"] == replicate
                ]

                experiment_list = replicate_subset["sample"].tolist()
                logging.debug("Samples: %s", experiment_list)

                remove_zeros(experiment_list, input_dir, output_dir, experiment_name)


if __name__ == "__main__":
    main()
