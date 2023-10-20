import csv
import logging
import pathlib
import os

import pandas as pd
import process_oligo_list
import process_variants
from Bio import SeqIO

# TODO: Expand indels to allow for greater than 3x!

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

experiment_name = snakemake.config["experiment"]


log_file = snakemake.log[0]
output_dir = "results/" + experiment_name + "/processed_counts/"
regenerate_variants = snakemake.params["regenerate_variants"]
remove_zeros = snakemake.params["remove_zeros"]


def process_experiment(
    experiment_list, regenerate_variants=False, remove_zeros=False
):
    """Process the experiment."""

    logging.debug("Processing files: %s", experiment_list)

    with open(os.path.join(snakemake.config["ref_dir"], snakemake.config["reference"]), "r") as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq

    if regenerate_variants:
        variant_list = process_oligo_list.designed_variants(
            snakemake.config["oligo_file"], str(ref_sequence), snakemake.config["offset"]
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
        df, other = process_variants.process_variants_file(csv_file, designed_df)

        # Write an Enrich2-readable output
        enrich_file = os.path.join(output_dir, experiment + ".tsv")
        process_variants.write_enrich_df(enrich_file, df)

        # Write the processed file to csv
        processed_file = os.path.join(output_dir, "counts", experiment + ".csv")
        p = pathlib.Path(processed_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(processed_file)

        # Write the rejected variants as well
        rejected_file = os.path.join(
            output_dir, "rejected", "rejected_" + experiment + ".csv"
        )
        p = pathlib.Path(rejected_file)
        p.parent.mkdir(parents=True, exist_ok=True)

        with p.open("w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(other)

    # Remove any variants with no observations before processing with Enrich2.

    if remove_zeros:
        enrich_file_list = [
            output_dir + s + ".tsv"
            for s in experiment_list
        ]
        zeros = process_variants.remove_zeros_enrich(enrich_file_list)

        unobserved_file = (
            output_dir + "/rejected/" + experiment_name + "_unobserved_variants.csv"
        )
        p = pathlib.Path(unobserved_file)
        p.parent.mkdir(parents=True, exist_ok=True)

        with p.open("w+") as f:
            for variant in zeros:
                f.write("%s\n" % variant)

        return

if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

logging.debug("Loading experiment file: %s", snakemake.config["experiment_file"])

with open(snakemake.config["experiment_file"]) as f:
    experiments = pd.read_csv(f, header=0).set_index(
        "sample", drop=False, verify_integrity=True
    )

for condition in experiments["condition"].unique():
    for replicate in experiments["replicate"].unique():
        experiment_list = experiments.loc[
            (experiments["condition"] == condition)
            & (experiments["replicate"] == replicate)
        ]["sample"].tolist()
        file_list = experiments.loc[
            (experiments["condition"] == condition)
            & (experiments["replicate"] == replicate)
        ]["file"].tolist()

        process_experiment(
            experiment_list, snakemake.params["regenerate_variants"], snakemake.params["remove_zeros"]
        )
