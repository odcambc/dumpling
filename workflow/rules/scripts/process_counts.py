import csv
import logging
import pathlib
import os

import pandas as pd
import process_oligo_list
import process_variants
from Bio import SeqIO

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

# Oligo processing files

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

log_file = snakemake.log[0]
output_dir = "results/" + snakemake.config["experiment"] + "/processed_counts/"

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

logging.debug("Loading experiment file: %s", snakemake.config["experiment_file"])


def process_experiment(
    experiment_list, regenerate_variants=False
):
    """Process the experiment."""

    logging.debug("Processing files: %s", experiment_list)

    with open(os.path.join(snakemake.config["ref_dir"], snakemake.config["reference"]), "r") as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq

    if regenerate_variants:
        offset = int(snakemake.config["orf"].split('-')[0]) - 4
        variant_list = process_oligo_list.designed_variants(
            snakemake.config["oligo_file"], str(ref_sequence), offset
        )
        process_oligo_list.write_designed_csv(
            snakemake.config["variants_file"], designed_variant_header, variant_list
        )

    with open(snakemake.config["variants_file"], "r") as f:
        variants_reader = csv.reader(f, delimiter=",")
        designed_variants_list = list(variants_reader)

    designed_df = pd.DataFrame.from_records(
        designed_variants_list[1:], columns=designed_variants_list[0]
    ).convert_dtypes()
    designed_df["pos"] = pd.to_numeric(designed_df["pos"])
    designed_df["count"] = pd.to_numeric(designed_df["count"])

    for experiment in experiment_list:
        csv_file = process_variants.read_gatk_csv(
            os.path.join(snakemake.params["gatk_dir"], experiment + ".variantCounts")
        )
        df, other, rejected_stats, accepted_stats, total_stats = process_variants.process_variants_file(csv_file, designed_df)

        # Write an Enrich2-readable output
        enrich_file = os.path.join(output_dir, experiment + ".tsv")
        process_variants.write_enrich_df(enrich_file, df)

        # Write the processed file to csv
        processed_file = os.path.join(output_dir, "counts", experiment + ".csv")
        p = pathlib.Path(processed_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(processed_file, index=False)

        # Write the rejected variants as well
        rejected_file = os.path.join(
            output_dir, "rejected", "rejected_" + experiment + ".csv"
        )
        p = pathlib.Path(rejected_file)
        p.parent.mkdir(parents=True, exist_ok=True)

        with p.open("w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(other)
        
        # Write the rejected stats file
        stats_file = os.path.join("stats", experiment_name, "processing", experiment + "_rejected_processing.tsv")
        process_variants.write_stats_file(stats_file, rejected_stats)

        # Write the accepted stats file
        stats_file = os.path.join("stats", experiment_name, "processing", experiment + "_accepted_processing.tsv")
        process_variants.write_stats_file(stats_file, accepted_stats)

        # Write the total stats file
        stats_file = os.path.join("stats", experiment_name, "processing", experiment + "_total_processing.tsv")
        process_variants.write_stats_file(stats_file, total_stats)

# Process the experiments according to their relationships
with open(snakemake.config["experiment_file"]) as f:
    experiments = pd.read_csv(f, header=0).dropna(how='all').set_index(
        "sample", drop=False, verify_integrity=True
    )

for condition in experiments["condition"].unique():
    logging.debug("Condition: %s", condition)
    tile_list = experiments.loc[
        (experiments["condition"] == condition)]["tile"]
    for tile in tile_list.unique():
        logging.debug("Tile: %s", tile)
        replicate_list = experiments.loc[
                (experiments["condition"] == condition)
                & (experiments["tile"] == tile)]["replicate"]
        for replicate in replicate_list.unique():
            logging.debug("Replicate: %s", replicate)
            experiment_name_list = experiments.loc[
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
            logging.debug(experiment_name_list)

            process_experiment(
                experiment_name_list, snakemake.params["regenerate_variants"]
            )
