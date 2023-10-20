import csv
import pathlib
import re

import regex
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


def name_to_hgvs(name):
    # syn case
    # insdel case
    return "p.(" + name + ")"


def designed_variants(oligo_csv, ref, offset):
    variant_list = []

    p = pathlib.Path(oligo_csv)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("r") as f:
        lines = csv.reader(f, delimiter=",")

        for line in lines:
            # Match substitutions
            variant = regex.search(
                r".*_DMS-[0-9]+_([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", line[0]
            )
            if variant:
                # Is this a synonymous mutation?
                if variant.group(1) == variant.group(3):
                    mutation_type = "S"
                else:
                    mutation_type = "M"

                codon = int(variant.groups(2)[1])

                pre_codon = ref[offset + 3 * codon - 9 : offset + 3 * codon]

                # Check and make sure it only splits once!
                # Add better failure handling
                # TODO: split is case-sensitive. We don't want that.

                ref_split = re.split(pre_codon, line[1], flags=re.IGNORECASE)
                if len(ref_split) == 2:
                    codon = ref_split[1][0:3]
                elif len(ref_split) == 1:
                    ref_split = ref_split = re.split(
                        pre_codon[2:], line[1], flags=re.IGNORECASE
                    )
                    try:
                        codon = ref_split[1][0:3]
                    except IndexError:
                        print(
                            "Error: ",
                            pre_codon,
                            codon,
                            ref_split,
                            variant,
                            mutation_type,
                        )

                name = (
                    seq1(variant.group(1)) + variant.group(2) + seq1(variant.group(3))
                )
                pos = int(variant.group(2))
                mutant = seq1(variant.group(3))
                length = 1

            # Deletions
            # The position of the deletion is set to be the *first* codon that is
            # deleted. Ex: F89_G91del would be a 3-codon deletion at position 89.

            variant = regex.search(r".*_delete-[0-9]+_([0-9]+)-([0-9]+)", line[0])
            if variant:
                mutation_type = "D"
                pos = int(variant.group(2))
                length = int(int(variant.group(1)) / 3)
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

            variant = regex.search(r".*_insert-[0-9]+_([a-zA-Z]+)-([0-9]+)", line[0])
            if variant:
                mutation_type = "I"
                pos = int(variant.group(2))
                length = int(len(variant.group(1)) / 3)
                start_codon = ref[offset + (3 * pos) : offset + 3 * (pos + 1)]
                start = Seq(start_codon).translate()
                end_codon = ref[offset + 3 * (pos + 1) : offset + 3 * (pos + 2)]
                end = Seq(end_codon).translate()
                if length == 1:
                    name = start + str(pos) + "_" + end + str(pos + 1) + "insG"
                elif length == 2:
                    name = start + str(pos) + "_" + end + str(pos + 1) + "insGS"
                elif length == 3:
                    name = start + str(pos) + "_" + end + str(pos + 1) + "insGSG"
                codon = ""
                mutant = "I_" + str(length)

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

    return variant_list


def write_designed_csv(file, header, variant_list):
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        csvwriter = csv.writer(f, delimiter=",")
        csvwriter.writerow(header)
        for row in variant_list:
            csvwriter.writerow(row)
