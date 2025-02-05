import csv
import logging
import pathlib
import os
import pandas as pd
from Bio import SeqIO

import process_oligo_list
import process_variants

from snakemake.script import snakemake  # Only used in main()


# Optional: Keep these lists as module-level constants if used by other functions
order: list[str] = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "D_1",
    "D_2",
    "D_3",
    "I_1",
    "I_2",
    "I_3",
    "X",
]

designed_variant_header: list[str] = [
    "count",
    "pos",
    "mutation_type",
    "name",
    "codon",
    "mutation",
    "length",
    "hgvs",
]


def process_experiment(
    experiment_list,
    ref_dir,
    reference_fasta,
    oligo_file,
    variants_file,
    orf_range,
    max_deletion_length,
    noprocess,
    gatk_dir,
    output_dir,
    regenerate_variants=False,
):
    """
    Process the experiment given a list of sample names (experiment_list).
    1. Load reference sequence from FASTA.
    2. Possibly regenerate variants file from oligo definitions.
    3. Read the designed variants file.
    4. For each sample in experiment_list:
       - Read GATK csv
       - Process variants
       - Write Enrich2-readable file, processed CSV, and stats.
    """

    logging.debug("Processing files: %s", experiment_list)

    # 1. Load reference sequence from FASTA
    ref_path = os.path.join(ref_dir, reference_fasta)
    with open(ref_path, "r") as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq

    # 2. Parse ORF range
    #    e.g. if orf_range == "100-500", orf_start=100, orf_end=500
    orf_start, orf_end = map(int, orf_range.split("-"))
    ref_AA_sequence = ref_sequence[orf_start - 1 : orf_end].translate()

    # 3. Possibly regenerate variants
    if regenerate_variants:
        offset = orf_start - 4
        variant_list = process_oligo_list.designed_variants(
            oligo_file, str(ref_sequence), offset
        )
        process_oligo_list.write_designed_csv(
            variants_file, designed_variant_header, variant_list
        )

    # 4. Read the designed variants file
    with open(variants_file, "r") as f:
        variants_reader = csv.reader(f, delimiter=",")
        designed_variants_list = list(variants_reader)

    designed_df = pd.DataFrame.from_records(
        designed_variants_list[1:], columns=designed_variants_list[0]
    ).convert_dtypes()
    designed_df["pos"] = pd.to_numeric(designed_df["pos"])
    designed_df["count"] = pd.to_numeric(designed_df["count"])

    # 5. Process each sample in the experiment list
    for experiment_name in experiment_list:
        # 5a. Read GATK CSV
        gatk_csv = os.path.join(gatk_dir, f"{experiment_name}.variantCounts")
        csv_file = process_variants.read_gatk_csv(gatk_csv)

        # 5b. Process variants (returns multiple dataframes & stats)
        df, other, rejected_stats, accepted_stats, total_stats = (
            process_variants.process_variants_file(
                csv_file, designed_df, ref_AA_sequence, max_deletion_length, noprocess
            )
        )

        # 5c. Write Enrich2-readable file
        enrich_file = os.path.join(
            output_dir, "enrich_format", f"{experiment_name}.tsv"
        )
        process_variants.write_enrich_df(enrich_file, df)

        # 5d. Write processed file to CSV
        processed_file = os.path.join(output_dir, f"{experiment_name}.csv")
        p = pathlib.Path(processed_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(processed_file, index=False)

        # 5e. Write rejected variants
        rejected_file = os.path.join(
            output_dir, "rejected", f"rejected_{experiment_name}.csv"
        )
        p = pathlib.Path(rejected_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(other)

        # 5f. Write rejected, accepted, and total stats files
        stats_dir = os.path.join("stats", snakemake.config["experiment"], "processing")
        pathlib.Path(stats_dir).mkdir(parents=True, exist_ok=True)

        # Rejected
        stats_file = os.path.join(
            stats_dir, f"{experiment_name}_rejected_processing.tsv"
        )
        process_variants.write_stats_file(stats_file, rejected_stats)

        # Accepted
        stats_file = os.path.join(
            stats_dir, f"{experiment_name}_accepted_processing.tsv"
        )
        process_variants.write_stats_file(stats_file, accepted_stats)

        # Total
        stats_file = os.path.join(stats_dir, f"{experiment_name}_total_processing.tsv")
        process_variants.write_stats_file(stats_file, total_stats)


def main():
    """
    Main function that orchestrates:
      - Logging setup
      - Reading config
      - Reading the experiment file
      - Looping over conditions, tiles, replicates
      - Calling process_experiment
    """

    # Set up logging
    log_file = snakemake.log[0]
    if log_file:
        logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG)

    # Read from Snakemake config
    experiment_name = snakemake.config["experiment"]
    experiment_file = snakemake.config["experiment_file"]
    ref_dir = snakemake.config["ref_dir"]
    reference_fasta = snakemake.config["reference"]
    orf_range = snakemake.config["orf"]  # e.g. "100-500"
    oligo_file = snakemake.config["oligo_file"]
    variants_file = snakemake.config["variants_file"]
    gatk_dir = snakemake.params["gatk_dir"]
    max_deletion_length = snakemake.config["max_deletion_length"]
    noprocess = snakemake.config["noprocess"]
    regenerate_variants = snakemake.params["regenerate_variants"]
    tiled = snakemake.config["tiled"]

    # Determine output directory
    output_dir = os.path.join("results", experiment_name, "processed_counts")

    logging.debug("Loading experiment file: %s", experiment_file)

    # Read experiments DataFrame
    with open(experiment_file) as f:
        experiments = (
            pd.read_csv(f, header=0)
            .dropna(how="all")
            .set_index("sample", drop=False, verify_integrity=True)
        )

    # Loop over conditions, tiles, replicates
    for condition in experiments["condition"].unique():
        logging.debug("Condition: %s", condition)

        if tiled:
            tile_list = experiments.loc[
                experiments["condition"] == condition, "tile"
            ].unique()
        else:
            # Treat as if there's only one tile when 'tile' doesn't exist or it's not tiled
            tile_list = [None]

        for tile in tile_list:
            logging.debug("Tile: %s", tile)
            if tile is not None:
                replicate_list = experiments.loc[
                    (experiments["condition"] == condition)
                    & (experiments["tile"] == tile),
                    "replicate",
                ]
            else:
                replicate_list = experiments.loc[
                    (experiments["condition"] == condition), "replicate"
                ]

            for replicate in replicate_list.unique():
                logging.debug("Replicate: %s", replicate)

                if tile is not None:
                    experiment_name_list = experiments.loc[
                        (experiments["condition"] == condition)
                        & (experiments["replicate"] == replicate)
                        & (experiments["tile"] == tile),
                        "sample",
                    ].tolist()
                else:
                    experiment_name_list = experiments.loc[
                        (experiments["condition"] == condition)
                        & (experiments["replicate"] == replicate),
                        "sample",
                    ].tolist()
                logging.debug("Samples: %s", experiment_name_list)

                # Call the processing function
                process_experiment(
                    experiment_list=experiment_name_list,
                    ref_dir=ref_dir,
                    reference_fasta=reference_fasta,
                    oligo_file=oligo_file,
                    variants_file=variants_file,
                    orf_range=orf_range,
                    max_deletion_length=max_deletion_length,
                    noprocess=noprocess,
                    gatk_dir=gatk_dir,
                    output_dir=output_dir,
                    regenerate_variants=regenerate_variants,
                )


if __name__ == "__main__":
    main()
