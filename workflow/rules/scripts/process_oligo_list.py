import csv
import pathlib
import re

import regex
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


PRE_SPAN = 15  # Number of bases in window to map


def name_to_hgvs(name):
    # syn case
    # insdel case
    return "p.(" + name + ")"


# TODO: add check if ref sequence is out of range.


def designed_variants(oligo_csv, ref, offset):
    variant_list = []

    p = pathlib.Path(oligo_csv)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("r") as f:
        lines = csv.reader(f, delimiter=",")

        for line in lines:
            # Match substitutions
            variant_sub = regex.search(
                r".*_DMS-[0-9]+_([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", line[0]
            )
            if variant_sub:
                # Is this a synonymous mutation?
                if variant_sub.group(1) == variant_sub.group(3):
                    mutation_type = "S"
                else:
                    mutation_type = "M"

                # The codon is assigned from processing the oligo itself.

                codon_n = int(variant_sub.groups(2)[1])

                pre_start = max(0, offset + (3 * codon_n) - PRE_SPAN)
                post_end = min(offset + (3 * codon_n) + PRE_SPAN, len(ref))

                pre_codon = ref[pre_start : offset + (3 * codon_n)]
                post_codon = ref[offset + (3 * codon_n) + 3 : post_end]

                pre_split = re.split(pre_codon, line[1], flags=re.IGNORECASE)

                if len(pre_split) == 2:
                    codon = pre_split[1][0:3]
                else:
                    post_split = re.split(post_codon, line[1], flags=re.IGNORECASE)
                    if len(post_split) == 2:
                        codon = post_split[0][-3:]
                    else:
                        print(
                            "Warning: could not find codon for",
                            line[0],
                            "Double check offset value and reference sequence.",
                        )

                name = (
                    seq1(variant_sub.group(1))
                    + variant_sub.group(2)
                    + seq1(variant_sub.group(3))
                )
                pos = int(variant_sub.group(2))
                mutant = seq1(variant_sub.group(3))
                length = 1

            # Deletions
            # The position of the deletion is set to be the *first* codon that is
            # deleted. Ex: F89_G91del would be a 3-codon deletion at position 89.
            # This assumes the length is given in *nucleotides*, not codons.

            variant_del = regex.search(r".*_delete-[0-9]+_([0-9]+)-([0-9]+)", line[0])

            if variant_del:
                mutation_type = "D"
                pos = int(variant_del.group(2))
                # check if length is divisible by 3
                # if not, print warning and set length to 0
                if int(variant_del.group(1)) % 3 != 0:
                    print("Warning: Deletion length is not divisible by 3")
                    length = 0
                else:
                    length = int(int(variant_del.group(1)) / 3)

                start_codon = ref[offset + (3 * pos) : offset + 3 * (pos + 1)]
                start = Seq(start_codon).translate()

                if length == 1:
                    name = start + str(pos) + "del"
                else:
                    end_codon = ref[
                        offset + 3 * (pos + length - 1) : offset + 3 * (pos + length)
                    ]
                    end = Seq(end_codon).translate()
                    name = start + str(pos) + "_" + end + str(pos + length - 1) + "del"
                codon = ""
                mutant = "D_" + str(length)

            # Insertions
            # The position of the insertion is defined as the position of the
            # first *inserted* residue. Ex: F47_V48insG is assigned position 48.

            variant_ins = regex.search(
                r".*_insert-[0-9]+_([a-zA-Z]+)-([0-9]+)", line[0]
            )
            if variant_ins:
                mutation_type = "I"
                pos = int(variant_ins.group(2))
                length = int(len(variant_ins.group(1)) / 3)
                start_codon = ref[offset + (3 * pos) : offset + 3 * (pos + 1)]
                start = Seq(start_codon).translate()
                end_codon = ref[offset + 3 * (pos + 1) : offset + 3 * (pos + 2)]
                end = Seq(end_codon).translate()

                codon = variant_ins.group(1)

                aa_string = ""

                for i in range(0, length):
                    aa_string = aa_string + str(
                        Seq(codon[(3 * i) : (3 * i) + 3]).translate()
                    )

                name = start + str(pos) + "_" + end + str(pos + 1) + "ins" + aa_string

                mutant = "I_" + str(length)

            try:
                variant_list = variant_list + [
                    [
                        0,
                        pos,
                        mutation_type,
                        name,
                        codon,
                        mutant,
                        length,
                        name_to_hgvs(name),
                    ]
                ]
            except UnboundLocalError:
                print("Error: no match found for", line[0])

    return variant_list


def write_designed_csv(file, header, variant_list):
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        csvwriter = csv.writer(f, delimiter=",")
        csvwriter.writerow(header)
        for row in variant_list:
            csvwriter.writerow(row)
