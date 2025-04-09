import pathlib
import csv
import regex
import re
import logging
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import seq1

import pandas as pd

from snakemake.script import snakemake

log_file = snakemake.log[0]

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

oligo_file = snakemake.input["oligo_file"]
variants_file = snakemake.output[0]

ref_file = snakemake.input["ref_fasta"]

# Define constants
PRE_SPAN = 15  # Number of bases in window to map


def name_to_hgvs(name):
    """
    Convert variant name to HGVS format.
    """
    return "p.(" + name + ")"


def get_sequence_segment(ref, start, end, is_circular=False):
    """
    Get a segment of the reference sequence, handling circular genomes if needed.

    Args:
        ref (str): Reference sequence
        start (int): Start position (0-based)
        end (int): End position (0-based, exclusive)
        is_circular (bool): Whether the reference is circular

    Returns:
        str: The requested segment of the reference sequence
    """
    ref_len = len(ref)

    if not is_circular:
        # For linear genomes, just return the segment directly
        return ref[max(0, start) : min(end, ref_len)]

    # For circular genomes, handle wrapping around the origin
    if start < 0 or end > ref_len:
        # Calculate positions that wrap around
        start_mod = start % ref_len
        end_mod = end % ref_len

        if start_mod < end_mod:
            # Simple case, no wrapping
            return ref[start_mod:end_mod]
        else:
            # Wrapping around the origin
            return ref[start_mod:] + ref[:end_mod]
    else:
        # No wrapping needed
        return ref[start:end]


def designed_variants(oligo_csv, ref, offset, is_circular=False):
    """
    Process a CSV file containing oligo sequences and extract mutation information.

    Args:
        oligo_csv (str): Path to CSV file containing oligo sequences
        ref (str): Reference DNA sequence
        offset (int): Offset to adjust codon positions
        is_circular (bool): Whether the reference is circular

    Returns:
        list: List of Variant named tuples containing mutation information
    """
    variant_list = []
    ref_len = len(ref)

    p = pathlib.Path(oligo_csv)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("r") as f:
        lines = csv.reader(f, delimiter=",")

        for i, line in enumerate(lines):
            if not line or len(line) < 2:  # Skip empty lines or insufficient columns
                logging.warning(f"Skipping line {i}: Insufficient data")
                continue

            # Skip header lines
            if line[0] == "name" or line[0] == "\ufeffname":
                logging.warning(f"Skipping header line {i}: {line[0]}")
                continue

            try:
                # Initialize variables to None to detect which pattern matched
                variant_data = {}

                # Match substitutions
                variant_sub = regex.search(
                    r".*_DMS-([0-9]+)_([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", line[0]
                )
                if variant_sub:
                    chunk = int(variant_sub.group(1))
                    # Is this a synonymous mutation?
                    if variant_sub.group(2) == variant_sub.group(4):
                        mutation_type = "S"  # Synonymous
                    else:
                        mutation_type = "M"  # Missense

                    codon_n = int(variant_sub.group(3))
                    codon_pos = offset + (
                        3 * (codon_n - 1)
                    )  # Adjust for 1-based position

                    pre_start = codon_pos - PRE_SPAN
                    post_end = codon_pos + 3 + PRE_SPAN

                    pre_codon = get_sequence_segment(
                        ref, pre_start, codon_pos + 3, is_circular
                    )
                    post_codon = get_sequence_segment(
                        ref, codon_pos + 6, post_end, is_circular
                    )

                    # Extract codon from oligo sequence
                    codon = extract_codon(
                        line[1],
                        pre_codon,
                        post_codon,
                        codon_pos,
                        ref,
                        offset,
                        codon_n,
                        is_circular,
                    )

                    name = (
                        seq1(variant_sub.group(2))
                        + variant_sub.group(3)
                        + seq1(variant_sub.group(4))
                    )
                    pos = int(variant_sub.group(3))
                    mutant = seq1(variant_sub.group(4))
                    length = 1

                    variant_data["count"] = 0
                    variant_data["pos"] = pos
                    variant_data["mutation_type"] = mutation_type
                    variant_data["name"] = name
                    variant_data["codon"] = codon
                    variant_data["mutant"] = mutant
                    variant_data["length"] = length
                    variant_data["hgvs"] = name_to_hgvs(name)
                    variant_data["chunk"] = chunk

                # Deletions
                # Updated regex to better handle deletion format variations
                variant_del = regex.search(
                    r".*_delete-([0-9]+)_([0-9]+)(?:-([0-9]+))?", line[0]
                )
                if variant_del and not variant_data:
                    mutation_type = "D"
                    deletion_chunk = int(variant_del.group(1))
                    deletion_length = int(variant_del.group(2))
                    pos = int(variant_del.group(3))

                    # Check if length is divisible by 3
                    if deletion_length % 3 != 0:
                        logging.warning(
                            f"Deletion length {deletion_length} is not divisible by 3 for {line[0]}"
                        )
                        codon_length = 0
                    else:
                        codon_length = deletion_length // 3

                    # Handle boundary cases
                    start_codon_pos = offset + (3 * (pos - 1))
                    start_codon = get_sequence_segment(
                        ref, start_codon_pos, start_codon_pos + 3, is_circular
                    )

                    if len(start_codon) == 3:
                        start = str(Seq(start_codon).translate())
                    else:
                        logging.warning(
                            f"Start codon incomplete ({len(start_codon)} bases) for {line[0]}"
                        )
                        start = "X"

                    if codon_length == 1:
                        name = start + str(pos) + "del"
                    else:
                        end_codon_pos = offset + 3 * ((pos + codon_length - 1) - 1)
                        end_codon = get_sequence_segment(
                            ref, end_codon_pos, end_codon_pos + 3, is_circular
                        )

                        if len(end_codon) == 3:
                            end = str(Seq(end_codon).translate())
                        else:
                            logging.warning(
                                f"End codon incomplete ({len(end_codon)} bases) for {line[0]}"
                            )
                            end = "X"

                        name = (
                            start
                            + str(pos)
                            + "_"
                            + end
                            + str(pos + codon_length - 1)
                            + "del"
                        )

                    codon = ""  # No specific codon for deletions
                    mutant = "D_" + str(codon_length)
                    length = codon_length

                    variant_data["count"] = 0
                    variant_data["pos"] = pos
                    variant_data["mutation_type"] = mutation_type
                    variant_data["name"] = name
                    variant_data["codon"] = codon
                    variant_data["mutant"] = mutant
                    variant_data["length"] = length
                    variant_data["hgvs"] = name_to_hgvs(name)
                    variant_data["chunk"] = deletion_chunk

                # Insertions
                variant_ins = regex.search(
                    r".*_insert-([0-9]+)_([a-zA-Z]+)-([0-9]+)", line[0]
                )
                if variant_ins and not variant_data:
                    mutation_type = "I"
                    insertion_chunk = int(variant_ins.group(1))
                    inserted_seq = variant_ins.group(2)
                    insertion_length = len(inserted_seq)
                    pos = int(variant_ins.group(3))

                    # Verify insertion length matches sequence
                    if insertion_length != len(inserted_seq):
                        logging.warning(
                            f"Insertion length {insertion_length} does not match sequence length {len(inserted_seq)} for {line[0]}"
                        )

                    length = len(inserted_seq) // 3
                    if len(inserted_seq) % 3 != 0:
                        logging.warning(
                            f"Insertion sequence length {len(inserted_seq)} is not divisible by 3 for {line[0]}"
                        )

                    # Handle boundary cases
                    start_codon_pos = offset + (3 * (pos - 1))
                    start_codon = get_sequence_segment(
                        ref, start_codon_pos, start_codon_pos + 3, is_circular
                    )

                    if len(start_codon) == 3:
                        start = str(Seq(start_codon).translate())
                    else:
                        logging.warning(
                            f"Start codon incomplete ({len(start_codon)} bases) for {line[0]}"
                        )
                        start = "X"

                    end_codon_pos = offset + (3 * pos)
                    end_codon = get_sequence_segment(
                        ref, end_codon_pos, end_codon_pos + 3, is_circular
                    )

                    if len(end_codon) == 3:
                        end = str(Seq(end_codon).translate())
                    else:
                        logging.warning(
                            f"End codon incomplete ({len(end_codon)} bases) for {line[0]}"
                        )
                        end = "X"

                    codon = inserted_seq
                    aa_string = ""

                    # Translate inserted sequence to amino acids
                    for i in range(length):
                        codon_seq = codon[(3 * i) : (3 * i) + 3]
                        if len(codon_seq) == 3:  # Ensure we have a complete codon
                            aa_string += str(Seq(codon_seq).translate())
                        else:
                            logging.warning(
                                f"Incomplete codon {codon_seq} in insertion for {line[0]}"
                            )
                            aa_string += "X"

                    name = (
                        start + str(pos) + "_" + end + str(pos + 1) + "ins" + aa_string
                    )
                    mutant = "I_" + str(length)

                    variant_data["count"] = 0
                    variant_data["pos"] = pos
                    variant_data["mutation_type"] = mutation_type
                    variant_data["name"] = name
                    variant_data["codon"] = codon
                    variant_data["mutant"] = mutant
                    variant_data["length"] = length
                    variant_data["hgvs"] = name_to_hgvs(name)
                    variant_data["chunk"] = insertion_chunk

                # If we found a match, add it to our list
                if variant_data:
                    variant_list.append(variant_data)
                else:
                    logging.warning(
                        f"No pattern matched for {line[0]}. This may be a header or unrecognized format."
                    )
                    logging.warning(f"Line: {line}")

            except Exception as e:
                logging.error(f"Error processing line {i}: {line[0]} - {str(e)}")
                # Log the full traceback for better debugging
                import traceback

                logging.error(traceback.format_exc())

    return variant_list


def extract_codon(
    oligo_sequence,
    pre_codon,
    post_codon,
    codon_pos,
    ref,
    offset,
    codon_n,
    is_circular=False,
):
    """
    Extract the codon from an oligo sequence using flanking sequences.

    Args:
        oligo_sequence (str): The oligo sequence
        pre_codon (str): Sequence before the codon
        post_codon (str): Sequence after the codon
        codon_pos (int): Position of the codon in the reference
        ref (str): Reference sequence
        offset (int): Offset value
        codon_n (int): Codon number
        is_circular (bool): Whether the reference is circular

    Returns:
        str: The extracted codon or empty string if not found
    """
    # Handle case where pre_codon or post_codon might be empty
    if not pre_codon and not post_codon:
        logging.warning(
            f"Both pre_codon and post_codon are empty for codon_n {codon_n}"
        )
        return ""

    if pre_codon:
        pre_split = re.split(pre_codon, oligo_sequence, flags=re.IGNORECASE)
        if len(pre_split) == 2:
            return pre_split[1][0:3]

    if post_codon:
        post_split = re.split(post_codon, oligo_sequence, flags=re.IGNORECASE)
        if len(post_split) == 2:
            return post_split[0][-3:]

    # Try alternative approach - direct lookup
    expected_codon = get_sequence_segment(ref, codon_pos, codon_pos + 3, is_circular)
    if expected_codon:
        expected_codon_pos = oligo_sequence.upper().find(expected_codon.upper())
        if expected_codon_pos >= 0:
            return oligo_sequence[expected_codon_pos : expected_codon_pos + 3]

    # If all attempts fail, log and return empty string
    logging.warning(f"Could not find codon. Check offset value and reference sequence.")
    logging.warning(
        f"codon_n: {codon_n}, pre_start: {codon_pos - PRE_SPAN}, post_end: {codon_pos + PRE_SPAN + 3}"
    )
    logging.warning(f"pre_codon: {pre_codon}, post_codon: {post_codon}")
    if pre_codon:
        logging.warning(
            f"pre_split: {re.split(pre_codon, oligo_sequence, flags=re.IGNORECASE)}"
        )
    if post_codon:
        logging.warning(
            f"post_split: {re.split(post_codon, oligo_sequence, flags=re.IGNORECASE)}"
        )
    return ""


def check_designed_df(df) -> bool:
    """
    Check the designed variants dataframe for any potential issues.

    Args:
        df (pd.DataFrame): Designed variants dataframe

    Returns:
        bool: True if the dataframe is valid, False otherwise
    """

    # Check for missing values
    missing_values = df.isnull().sum().sum()
    if missing_values > 0:
        logging.error(
            f"Found {missing_values} missing values in the designed variants dataframe"
        )
        return False

    return True


def write_designed_csv(file, header, variant_list):
    """
    Write variants to a CSV file.

    Args:
        file (str): Output file path
        header (list): Column headers
        variant_list (list): List of variant dicts
    """
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        csvwriter = csv.writer(f, delimiter=",")
        csvwriter.writerow(header)
        for row in variant_list:
            csvwriter.writerow(row)


# Main execution

orf_range = snakemake.config[
    "orf"
]  # e.g. "100-500" or "500-100" for circular crossing origin

orf_parts = orf_range.split("-")
orf_start = int(orf_parts[0])
orf_end = int(orf_parts[1])

# Determine if the reference is circular based on the ORF range
is_circular = orf_start > orf_end
if is_circular:
    logging.info("Detected circular genome (ORF crosses origin)")
    # For circular genomes with reversed ranges, we'll adjust the offset differently
    offset = orf_start - 4  # Maintain the same offset calculation for consistency
else:
    offset = orf_start - 4  # Adjust offset based on ORF start

# Read in the reference sequence
with open(ref_file, "r") as f:
    ref_list = list(SeqIO.parse(f, "fasta"))
    ref_sequence = ref_list[0].seq

# For a circular genome, if ORF crosses the origin, we need to handle the ORF extraction differently
if is_circular:
    # Extract ORF that wraps around the origin
    # First part: from orf_start to the end of sequence
    first_part = ref_sequence[orf_start - 1 :]
    # Second part: from the beginning to orf_end
    second_part = ref_sequence[:orf_end]
    # Join them to get the complete ORF
    ref_AA_sequence = (first_part + second_part).translate()
    logging.info(
        f"Circular ORF: joined sequence from positions {orf_start-1}:{len(ref_sequence)} and 0:{orf_end}"
    )
else:
    # Linear case - just extract the ORF directly
    ref_AA_sequence = ref_sequence[orf_start - 1 : orf_end].translate()

logging.info(
    f"ORF range: {orf_range}, offset: {offset}, ORF start: {orf_start}, ORF end: {orf_end}, is_circular: {is_circular}"
)

logging.info(f"Reference sequence length: {len(ref_sequence)}")
logging.info(f"Reference AA sequence length: {len(ref_AA_sequence)}")

# Generate the variants: returns list of dicts
variants = designed_variants(oligo_file, str(ref_sequence), offset, is_circular)
logging.info(variants[0:20])
logging.info(f"Generated {len(variants)} variants.")

# Convert the list of dicts to a DataFrame
variants_df = pd.DataFrame(variants)
# Make sure types are strings
variants_df["name"] = variants_df["name"].astype(str)
variants_df["codon"] = variants_df["codon"].astype(str)
variants_df["mutant"] = variants_df["mutant"].astype(str)
variants_df["hgvs"] = variants_df["hgvs"].astype(str)


# Check for any issues in the designed variants
if not check_designed_df(variants_df):
    logging.error("Error in designed variants. Check log for details.")
    raise Exception("Error in designed variants. Check log for details.")

# Check for duplicate rows
duplicate_rows = variants_df.duplicated().sum()

logging.info(variants_df.nunique())
logging.info(
    f"Found {duplicate_rows} completely duplicate rows in the designed variants dataframe."
)

# Deduplicate identical rows
if duplicate_rows > 0:
    logging.warning(
        f"Found {duplicate_rows} duplicate rows in the designed variants dataframe. Dropping duplicates."
    )
    variants_df = variants_df.drop_duplicates()

# If there are duplicated names without identical positions, etc, something else might be wrong.
duplicated_names = variants_df.drop_duplicates()["name"].duplicated().sum()
if duplicated_names > 0:
    logging.warning(
        f"Found {duplicated_names} duplicated variant names in the designed variants dataframe. Attempting to deduplicate."
    )

    # Try to merge duplicates. If two rows have the same name but only differ in chunk values, we can combine and set the chunk to the first.
    # In this case, the codon might also be different, if the variant was generated in two different chunks!
    # We will keep the first chunk and codon, and drop the others.
    variants_df = (
        variants_df.groupby(
            [
                "count",
                "pos",
                "mutation_type",
                "name",
                "mutant",
                "length",
                "hgvs",
            ],
            as_index=False,
        )
        .agg({"chunk": "first", "codon": "first"})
        .reset_index()
    )

    # Check again for duplicates
    duplicated_names = variants_df["name"].duplicated().sum()
    if duplicated_names > 0:
        # If we still have duplicates, it means they are not trivially different.
        variants_df.to_csv("duped.csv", index=False)
        logging.error(
            f"Found {duplicated_names} duplicated variant names with non-identical values in the designed variants dataframe. Check for errors."
        )
        raise Exception(
            "Found duplicated variant names with non-identical values. Check for errors."
        )

# Write the variants to a file
logging.info("Regenerated variants.")

variants_df.to_csv(variants_file, index=False)

logging.info(f"New designed variants file written to {variants_file}.")
