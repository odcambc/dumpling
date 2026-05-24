import csv
import logging
import pathlib
import os
import pandas as pd
from Bio import SeqIO

import process_variants
from script_utils import run_script, set_index_with_unique_check


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
    "chunk",
]


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


def write_stats_output(stats_file, stats):
    """
    Write stats to a file.
    """

    with open(stats_file, "w") as f:
        f.write("Mutation\tCount\n")
        for mutation in order:
            f.write(f"{mutation}\t{stats.get(mutation, 0)}\n")


def write_counts_output(counts_file, counts):
    """
    Write counts to a file.
    """

    with open(counts_file, "w") as f:
        f.write("Mutation\tCount\n")
        for mutation in order:
            f.write(f"{mutation}\t{counts.get(mutation, 0)}\n")


def process_experiment(
    sample_list,
    experiment_name,
    ref_AA_sequence,
    designed_df,
    max_deletion_length,
    noprocess,
    gatk_dir,
    output_dir,
):
    """
    Process one replicate-group's worth of samples.

    The reference AA sequence and designed-variants DataFrame are passed in
    by `_run` rather than loaded here. Both are identical across every
    (condition × tile × replicate) call, so the FASTA parse + ORF translate
    and the variants CSV read happen once per pipeline invocation instead
    of once per replicate group (audit item M7).

    For each sample in `sample_list`:
      - Read GATK csv
      - Process variants
      - Write Enrich2-readable file, processed CSV, and stats.
    """

    logging.debug("Processing files: %s", sample_list)

    # Process each sample in the experiment list
    for sample_name in sample_list:
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

        # 5d. Write processed file to CSV
        processed_file = os.path.join(output_dir, f"{sample_name}.csv")
        p = pathlib.Path(processed_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        logging.debug(
            f"Found {len(filtered_df)} accepted variants out of {len(designed_df)} possible variants"
        )
        filtered_df.to_csv(processed_file, index=False)

        # 5e. Write rejected variants
        rejected_file = os.path.join(
            output_dir, "rejected", f"rejected_{sample_name}.csv"
        )
        p = pathlib.Path(rejected_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(rejected_df)

        # 5f. Write rejected, accepted, and total stats files
        stats_dir = os.path.join("stats", experiment_name, "processing")
        pathlib.Path(stats_dir).mkdir(parents=True, exist_ok=True)

        # Rejected
        stats_file = os.path.join(stats_dir, f"{sample_name}_rejected_processing.tsv")
        process_variants.write_stats_file(stats_file, rejected_stats)

        # Accepted
        stats_file = os.path.join(stats_dir, f"{sample_name}_accepted_processing.tsv")
        process_variants.write_stats_file(stats_file, accepted_stats)

        # Total
        stats_file = os.path.join(stats_dir, f"{sample_name}_total_processing.tsv")
        process_variants.write_stats_file(stats_file, total_stats)


def main():
    from snakemake.script import snakemake

    run_script(snakemake, _run)


def _run(snakemake):
    # Read from Snakemake config
    experiment_name = snakemake.config["experiment"]
    experiment_file = snakemake.config["experiment_file"]
    ref_dir = snakemake.config["ref_dir"]
    reference_fasta = snakemake.config["reference"]
    orf_range = snakemake.config["orf"]  # e.g. "100-500"
    designed_variants_file = snakemake.config["variants_file"]
    gatk_dir = snakemake.params["gatk_dir"]
    noprocess = snakemake.config["noprocess"]
    max_deletion_length = snakemake.config["max_deletion_length"]
    tiled = snakemake.config["tiled"]

    # Determine output directory
    output_dir = os.path.join("results", experiment_name, "processed_counts")

    logging.debug("Loading experiment file: %s", experiment_file)

    # Read experiments DataFrame
    with open(experiment_file, encoding="utf-8-sig") as f:
        samples = set_index_with_unique_check(
            pd.read_csv(f, header=0).dropna(how="all"),
            "sample",
            drop=False,
        )

    # Parse reference FASTA + translate ORF region once, then reuse across
    # every (condition × tile × replicate) call below. Audit item M7:
    # these were previously redone per-call inside process_experiment.
    ref_path = os.path.join(ref_dir, reference_fasta)
    with open(ref_path, "r") as f:
        ref_list = list(SeqIO.parse(f, "fasta"))
        ref_sequence = ref_list[0].seq
    orf_start, orf_end = map(int, orf_range.split("-"))
    ref_AA_sequence = ref_sequence[orf_start - 1 : orf_end].translate()

    # Read designed variants file once. When noprocess=True the GATK output
    # is used as-is (no filtering); pass an empty frame through and skip
    # the file existence check.
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

    # Loop over conditions, tiles, replicates
    for condition in samples["condition"].unique():
        logging.debug("Condition: %s", condition)

        if tiled:
            tile_list = samples.loc[samples["condition"] == condition, "tile"].unique()
        else:
            # Treat as if there's only one tile when 'tile' doesn't exist or it's not tiled
            tile_list = [None]

        for tile in tile_list:
            logging.debug("Tile: %s", tile)
            if tile is not None:
                replicate_list = samples.loc[
                    (samples["condition"] == condition) & (samples["tile"] == tile),
                    "replicate",
                ]
            else:
                replicate_list = samples.loc[
                    (samples["condition"] == condition), "replicate"
                ]

            for replicate in replicate_list.unique():
                logging.debug("Replicate: %s", replicate)

                if tile is not None:
                    sample_name_list = samples.loc[
                        (samples["condition"] == condition)
                        & (samples["replicate"] == replicate)
                        & (samples["tile"] == tile),
                        "sample",
                    ].tolist()
                else:
                    sample_name_list = samples.loc[
                        (samples["condition"] == condition)
                        & (samples["replicate"] == replicate),
                        "sample",
                    ].tolist()
                logging.debug("Samples: %s", sample_name_list)

                # Call the processing function
                process_experiment(
                    sample_list=sample_name_list,
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
