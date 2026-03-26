import logging
from pathlib import Path

import pandas as pd

from script_utils import run_script, set_index_with_unique_check


def remove_zeros_enrich(enrich_file_list, output_dir):
    """
    Combine Enrich2 TSV files, remove unobserved variants across all experiments,
    and rewrite the filtered TSV files.

    :param enrich_file_list: List of input file paths in Enrich2 format
    :param output_dir: Directory to place filtered output
    :return: List of unobserved variant IDs
    """
    # Build a combined DataFrame from all inputs

    df_list = []
    for enrich_file in enrich_file_list:
        file_stem = Path(enrich_file).stem  # e.g. "sample1"
        # Use "Int64" so we can safely store 0/NA as integers
        df = pd.read_csv(
            enrich_file,
            sep="\t",
            names=["hgvs", file_stem],
            header=0,
            dtype={"hgvs": str, file_stem: "Int64"},
        ).set_index("hgvs")
        df_list.append(df)

    combined_enrich_df = pd.concat(df_list, axis=1, join="outer")
    # NA or 0 means unobserved. We define a boolean mask:
    is_null_or_zero = combined_enrich_df.isna() | (combined_enrich_df == 0)

    # Variants that are unobserved in ALL files
    unobserved_mask = is_null_or_zero.all(axis=1)
    unobserved_variants = unobserved_mask[unobserved_mask].index.to_list()

    # Keep only rows that have at least one non-zero, non-NA value
    combined_enrich_df = combined_enrich_df[~unobserved_mask]

    # Rewrite each file using only the "surviving" variants
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for enrich_file in enrich_file_list:
        file_stem = Path(enrich_file).stem
        file_name = Path(enrich_file).name
        output_path = output_dir / file_name

        # We write just the single column for this file
        # The index is 'hgvs', and the column is the newly filtered data
        # named after the file_stem
        subset_df = combined_enrich_df[[file_stem]].fillna(0).astype("Int64")

        # Rename the column back to "count" so it remains in Enrich2 format
        subset_df.columns = ["count"]

        # Write out
        subset_df.to_csv(
            output_path,
            header=True,
            index=True,
            sep="\t",
        )
        logging.debug("Writing zero-filtered file %s -> %s", enrich_file, output_path)

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


def _run(snakemake):
    """
    Main function to perform unobserved variant removal for each condition/tile/replicate grouping.
    - Reads config from Snakemake
    - Checks if 'tile' column exists for optional tiling
    - Calls remove_zeros() for each replicate group
    """
    logging.debug("Performing unobserved variant removal.")

    experiment_name = snakemake.config["experiment"]
    experiment_file = snakemake.config["experiment_file"]
    input_dir = f"results/{experiment_name}/processed_counts/enrich_format/"
    output_dir = f"results/{experiment_name}/processed_counts/removed_zeros/"

    # Read experiments
    experiments = set_index_with_unique_check(
        pd.read_csv(experiment_file, header=0).dropna(how="all"),
        "sample",
        drop=False,
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


def main():
    from snakemake.script import snakemake

    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
