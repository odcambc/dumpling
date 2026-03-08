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


def setup_logging(log_file):
    if log_file:
        logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG)


def name_to_hgvs(name):
    return "p.(" + name + ")"


def get_sequence_segment(ref, sequence, start, end, is_circular=False):
    ref_len = len(ref)
    if not is_circular:
        return ref[max(0, start) : min(end, ref_len)]
    if start < 0 or end > ref_len:
        start_mod, end_mod = start % ref_len, end % ref_len
        return (
            ref[start_mod:] + ref[:end_mod]
            if start_mod > end_mod
            else ref[start_mod:end_mod]
        )
    return ref[start:end]


def parse_substitution(name, sequence, ref, offset, is_circular):
    match = regex.search(r".*_DMS-([0-9]+)_([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", name)
    if not match:
        return None
    chunk, original, position, mutant = match.groups()
    mutation_type = "S" if original == mutant else "M"
    codon_n = int(position)
    codon_pos = offset + (3 * (codon_n - 1))
    pre_codon = get_sequence_segment(ref, codon_pos - 15, codon_pos + 3, is_circular)
    post_codon = get_sequence_segment(ref, codon_pos + 6, codon_pos + 18, is_circular)
    name = seq1(original) + position + seq1(mutant)
    return [
        0,
        int(position),
        mutation_type,
        name,
        "",
        mutant,
        1,
        name_to_hgvs(name),
        int(chunk),
    ]


def parse_deletion(name, sequence, ref, offset, is_circular):
    match = regex.search(r".*_delete-([0-9]+)_([0-9]+)-([0-9]+)?", name)
    if not match:
        return None
    chunk, deletion_length, pos = map(int, match.groups())
    if deletion_length % 3 != 0:
        logging.warning(f"Deletion length {deletion_length} is not divisible by 3")
        return None
    codon_length = deletion_length // 3
    start_codon_pos = offset + (3 * (pos - 1))
    start = str(
        Seq(
            get_sequence_segment(ref, start_codon_pos, start_codon_pos + 3, is_circular)
        ).translate()
    )
    name = (
        f"{start}{pos}del"
        if codon_length == 1
        else f"{start}{pos}_X{pos + codon_length - 1}del"
    )
    return [
        0,
        pos,
        "D",
        name,
        "",
        f"D_{codon_length}",
        codon_length,
        name_to_hgvs(name),
        chunk,
    ]


def parse_insertion(name, sequence, ref, offset, is_circular):
    match = regex.search(r".*_insert-([0-9]+)_([a-zA-Z]+)-([0-9]+)", name)
    if not match:
        return None
    chunk, inserted_seq, pos = match.groups()
    length = len(inserted_seq) // 3
    start_codon_pos = offset + (3 * (int(pos) - 1))
    start = Seq(
        get_sequence_segment(ref, start_codon_pos, start_codon_pos + 3, is_circular)
    ).translate()
    aa_string = "".join(
        str(Seq(inserted_seq[i : i + 3]).translate())
        for i in range(0, len(inserted_seq), 3)
    )
    name = f"{start}{pos}_X{int(pos) + 1}ins{aa_string}"
    return [
        0,
        pos,
        "I",
        name,
        inserted_seq,
        f"I_{length}",
        length,
        name_to_hgvs(name),
        int(chunk),
    ]


def designed_variants(oligo_csv, ref, offset, is_circular):
    variant_list = []
    with open(oligo_csv, "r") as f:
        lines = csv.reader(f, delimiter=",")
        for i, line in enumerate(lines):
            if not line or len(line) < 2:
                logging.warning(f"Skipping line {i}: Insufficient data")
                continue
            name, sequence = line[:2]
            for parser in (parse_substitution, parse_deletion, parse_insertion):
                variant = parser(name, sequence, ref, offset, is_circular)
                if variant:
                    variant_list.append(variant)
                    break
            else:
                logging.warning(f"No pattern matched for {name}. Skipping.")
    return variant_list


def write_designed_csv(file, header, variant_list):
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w+") as f:
        csv.writer(f, delimiter=",").writerow(header)
        csv.writer(f, delimiter=",").writerows(variant_list)


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
        f"Circular ORF: joined sequence from positions {orf_start - 1}:{len(ref_sequence)} and 0:{orf_end}"
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
]

write_designed_csv(variants_file, header, variants)
