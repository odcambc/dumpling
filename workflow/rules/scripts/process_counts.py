import csv
import logging
import pathlib
import os
import pandas as pd
from Bio import SeqIO

import process_variants
from script_utils import run_script, translate_orf


def process_gatk_file(gatk_output_file, designed_df, ref_AA_sequence, max_deletion_length, noprocess):
    """
    Process a GATK CSV file.
    * Read the GATK CSV file.
    * Process the variants.
    * Write the Enrich2-readable file, processed CSV, and stats.
    """

    # Read GATK CSV
    gatk_data = process_variants.read_gatk_csv(gatk_output_file)

    # Process variants (returns multiple dataframes & stats)
    filtered_df, rejected_df, rejected_stats, accepted_stats, total_stats = (
        process_variants.process_variants_file(
            gatk_data, designed_df, ref_AA_sequence, max_deletion_length, noprocess
        )
    )

    return filtered_df, rejected_df, rejected_stats, accepted_stats, total_stats


def process_sample(
    sample_name,
    experiment_name,
    ref_AA_sequence,
    designed_df,
    max_deletion_length,
    noprocess,
    gatk_dir,
    output_dir,
):
    """
    Process one sample's GATK output: filter variants, write Enrich2 + CSV + stats.

    Snakemake invokes this once per sample via the `process_sample` rule's
    `sample_prefix` wildcard, so each call is its own Python process and can
    run in parallel with siblings on a cluster or multi-core local host.
    """

    logging.debug("Processing sample: %s", sample_name)

    gatk_output_file = os.path.join(gatk_dir, f"{sample_name}.variantCounts")
    logging.debug("Reading GATK file: %s", gatk_output_file)

    filtered_df, rejected_df, rejected_stats, accepted_stats, total_stats = (
        process_gatk_file(gatk_output_file, designed_df, ref_AA_sequence, max_deletion_length, noprocess)
    )
    logging.debug("Finished processing GATK file")

    # Write Enrich2-readable file
    enrich_file = os.path.join(output_dir, "enrich_format", f"{sample_name}.tsv")
    logging.debug(f"Writing Enrich2 file: {enrich_file}")

    process_variants.write_enrich_df(enrich_file, filtered_df, noprocess)

    # Write processed file to CSV
    processed_file = os.path.join(output_dir, f"{sample_name}.csv")
    p = pathlib.Path(processed_file)
    p.parent.mkdir(parents=True, exist_ok=True)
    logging.debug(
        f"Found {len(filtered_df)} accepted variants out of {len(designed_df)} possible variants"
    )
    filtered_df.to_csv(processed_file, index=False)

    # Write rejected variants
    rejected_file = os.path.join(
        output_dir, "rejected", f"rejected_{sample_name}.csv"
    )
    p = pathlib.Path(rejected_file)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerows(rejected_df)

    # Write rejected, accepted, and total stats files
    stats_dir = os.path.join("stats", experiment_name, "processing")
    pathlib.Path(stats_dir).mkdir(parents=True, exist_ok=True)

    stats_file = os.path.join(stats_dir, f"{sample_name}_rejected_processing.tsv")
    process_variants.write_stats_file(stats_file, rejected_stats)

    stats_file = os.path.join(stats_dir, f"{sample_name}_accepted_processing.tsv")
    process_variants.write_stats_file(stats_file, accepted_stats)

    stats_file = os.path.join(stats_dir, f"{sample_name}_total_processing.tsv")
    process_variants.write_stats_file(stats_file, total_stats)


def main():
    run_script(snakemake, _run)


def _run(snakemake):
    # Read from Snakemake config and wildcards
    experiment_name = snakemake.config["experiment"]
    ref_dir = snakemake.config["ref_dir"]
    reference_fasta = snakemake.config["reference"]
    orf_range = snakemake.config["orf"]  # e.g. "100-500"
    designed_variants_file = snakemake.config["variants_file"]
    gatk_dir = snakemake.params["gatk_dir"]
    noprocess = snakemake.config["noprocess"]
    max_deletion_length = snakemake.config["max_deletion_length"]
    sample_name = snakemake.wildcards.sample_prefix

    # Determine output directory
    output_dir = os.path.join("results", experiment_name, "processed_counts")

    # Parse reference FASTA + translate ORF region.
    ref_path = os.path.join(ref_dir, reference_fasta)
    with open(ref_path, "r") as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq
    ref_AA_sequence = translate_orf(ref_sequence, orf_range)

    # Read designed variants file. When noprocess=True the GATK output is used
    # as-is (no filtering); pass an empty frame through and skip the file
    # existence check.
    if noprocess:
        logging.info("noprocess=True: skipping designed variants file load.")
        designed_df = pd.DataFrame()
    else:
        logging.info("Loading designed variants file: %s", designed_variants_file)
        if not os.path.exists(designed_variants_file):
            raise FileNotFoundError(
                f"Variants file not found: '{designed_variants_file}'. "
                "Check 'variants_file' path in config YAML."
            )
        designed_df = pd.read_csv(designed_variants_file, encoding="utf-8-sig")
        logging.info("Designed variants length: %d", len(designed_df))

    process_sample(
        sample_name=sample_name,
        experiment_name=experiment_name,
        ref_AA_sequence=ref_AA_sequence,
        designed_df=designed_df,
        max_deletion_length=max_deletion_length,
        noprocess=noprocess,
        gatk_dir=gatk_dir,
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
