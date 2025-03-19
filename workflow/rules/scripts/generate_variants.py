import pathlib
import csv
import regex
import re
import logging
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import seq1
from collections import namedtuple

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

# Define a named tuple for variant data
Variant = namedtuple(
    "Variant",
    [
        "index",
        "position",
        "mutation_type",
        "name",
        "codon",
        "mutant",
        "length",
        "hgvs",
        "chunk",
    ],
)


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
        list: List of Variant namedtuples containing mutation information
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

            try:
                # Initialize variables to None to detect which pattern matched
                variant_data = None

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

                    variant_data = [
                        0,
                        pos,
                        mutation_type,
                        name,
                        codon,
                        mutant,
                        length,
                        name_to_hgvs(name),
                        chunk,
                    ]

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

                    variant_data = [
                        0,
                        pos,
                        mutation_type,
                        name,
                        codon,
                        mutant,
                        length,
                        name_to_hgvs(name),
                        deletion_chunk,
                    ]

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

                    variant_data = [
                        0,
                        pos,
                        mutation_type,
                        name,
                        codon,
                        mutant,
                        length,
                        name_to_hgvs(name),
                        insertion_chunk,
                    ]

                # If we found a match, add it to our list
                if variant_data:
                    variant_list.append(Variant(*variant_data))
                else:
                    logging.warning(
                        f"No pattern matched for {line[0]}. This may be a header or unrecognized format."
                    )

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


def write_designed_csv(file, header, variant_list):
    """
    Write variants to a CSV file.

    Args:
        file (str): Output file path
        header (list): Column headers
        variant_list (list): List of Variant namedtuples
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

# Generate the variants
variants = designed_variants(oligo_file, str(ref_sequence), offset, is_circular)

# Write the variants to a file
header = [
    "index",
    "position",
    "mutation_type",
    "name",
    "codon",
    "mutant",
    "length",
    "hgvs",
    "chunk",
]

write_designed_csv(variants_file, header, variants)
