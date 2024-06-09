import csv
import pathlib
import logging

import pandas as pd
import regex

# TODO: create stats file and return

aa_3to1_dict = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Stp": "X",
}
aa_1to3_dict = {v: k for k, v in aa_3to1_dict.items()}

variantCounts_colnames = [
    "counts",
    "coverage",
    "mean_length",
    "length_NT",
    "NT",
    "length_codon",
    "codon",
    "AA",
    "mutations",
]

def name_to_hgvs(name):
    # syn case
    # insdel case
    return "p.(" + name + ")"


def read_gatk_csv(file):
    # GATK output unfortunately may not be full-width, which will confuse pandas during
    # reading in. We'll read each line and make sure it's full-length before passing it
    # to pandas.
    logging.info("Loading GATK output file: %s", file)

    gatk_list = []

    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("r") as f:
        lines = csv.reader(f, delimiter="\t")
        for line in lines:
            while len(line) < 9:
                line = line + [""]

            gatk_list.append(line)

    return gatk_list


# Process variant types
# The following functions process specific types of variants obtained from
# GATK ASM. Each accepts input a line from ASM output and parses it, and returns
# a dict containing variant info:
#           int count - number of observations of variant
#           int pos - position of variant (in codons)
#           str mutation_type - mutation class: S = synonymous
#                                      M = substitution
#                                      I = insertion
#                                      D = deletion
#           str name - mutation name, in 1 letter amino acids
#           str codon - for substitutions, the designed codon. For insertions, the added codons.
#                       empty for others
#           str mutation - for substitutions, the variant residue. for
#                           indels, the specific variant, in the form (I/D)_(1/2/N),
#                           as appropriate
#           int length - length of change, in codons
#           int rejected - whether or not the variant should be rejected


def process_insertion(line):
    # This contains the logic for parsing insertion variants of arbitrary length.
    # We define the position of an insertion as the position of the first new codon
    # in the new frame.

    variant_dict = {}

    count = int(line[0])
    length = int(line[5])
    mutation_type = "I"
    variant = line[8]
    mutation = "I_" + str(length)
    codons = line[6]
    pos = int(codons.split(":")[0])
    rejected = 0

    codon = ""

    if length > 1:
        for c in codons.split(", "):
            codon = codon + c.split(">")[1]
    else:
        codon = codon + codons.split(">")[1]

    variant_dict["counts"] = count
    variant_dict["pos"] = pos
    variant_dict["mutation_type"] = mutation_type
    variant_dict["name"] = variant
    variant_dict["codon"] = codon
    variant_dict["mutation"] = mutation
    variant_dict["length"] = length
    variant_dict["rejected"] = rejected

    return variant_dict


def process_deletion(line):
    # This contains the logic for parsing deletion variants of arbitrary length.
    # We define the position of deletion as the position of the first deleted codon.

    variant_dict = {}
    codon = line[6]
    length_codon = int(line[5])
    variant = line[8]
    rejected = 0

    if "_" in variant:
        start_codon = int(variant.split("_")[0][1:])
        end_codon = int(variant.split("_")[1][1:-3])

        # It may be the case that deletions are coupled with synonymous variants.
        # Reject these.
        if length_codon != (end_codon - start_codon + 1):
            rejected = 1
    elif length_codon > 1:
        rejected = 1

    variant_dict["counts"] = int(line[0])
    variant_dict["pos"] = int(codon.split(":")[0])
    variant_dict["mutation_type"] = "D"
    variant_dict["name"] = variant
    variant_dict["codon"] = ""
    variant_dict["mutation"] = "D_" + str(length_codon)
    variant_dict["length"] = length_codon
    variant_dict["rejected"] = rejected

    return variant_dict


def process_insdel(line, variant_names_array, ref_AA_sequence, max_deletion_length):
    # This contains the logic for parsing insdel variants and checking for an edge
    # case where in-frame deletions are being called as insdels.
    # Deletions of the form XYZA->X--- (just deletion) may be called as
    # XYZA->---X (deletion + mutation)
    # Note: bbmap defaults to left alignment, as is standard for NGS reads.
    # This means we only need to consider incorrectly left-aligned cases.

    # The logic here requires passing an array of the designed variant names to
    # search for the proper canonical name. It also requires returning an additional
    # "rejected" value in the dict.

    variant_dict = {}
    codon = line[6]
    length_codon = int(line[5])
    variant = line[8]
    count = int(line[0])
    rejected = 1
    mutation_type = "Z"
    mutation = "Z"
    pos = -1
    name = ""

    # Matches things like G37_I38insdelG or H296_V297insdelGS
    # Groups are: (1) start AA, (2) start pos
    # (3) end AA, (4) end pos, (5) inserted AAs
    # length is end pos - start pos (length of "deleted" region)
    # insdel_length is number of "inserted" AAs

    insdel_re = regex.search(
        r"([a-zA-Z])([0-9]+)_([a-zA-Z])([0-9]+)insdel([a-zA-Z]+)", variant
    )

    start_aa = insdel_re.group(1)
    end_aa = insdel_re.group(3)
    start_pos = int(insdel_re.group(2))
    end_pos = int(insdel_re.group(4))
    insdel_aas = insdel_re.group(5)

    deletion_length = end_pos - start_pos + 1  # Total number of AAs deleted
    insertion_length = len(insdel_aas)  # Total number of AAs inserted

    insdel_length = (
        deletion_length - insertion_length
    )  # Effective number of deleted AAs

    if (
        start_aa == insdel_aas  # start and inserted AA are the same
        or end_aa == insdel_aas  # end and inserted AA are the same
        or start_aa == end_aa  # start and end AA are the same
    ):
        if (
            insdel_length <= max_deletion_length
        ):  # the deletion should be an expected length
            rejected = 0
            pos = start_pos + insertion_length  # position of first changed AA
            mutation_type = "D"
            mutation = "D_" + str(insertion_length)

            if insdel_length == 1:
                name = end_aa + str(pos) + "del"
            elif insdel_length > 1:
                name = (
                    ref_AA_sequence[pos - 1]
                    + str(pos)
                    + "_"
                    + ref_AA_sequence[pos + insdel_length - 2]
                    + str(pos + insdel_length - 1)
                    + "del"
                )
            else:
                name = ""
                logging.warning(
                    "Error in insdel length calculation: insdel_length %i",
                    insdel_length,
                )
                rejected = 1

            if name in variant_names_array:
                rejected = 0
            else:
                rejected = 1

    variant_dict["counts"] = count
    variant_dict["pos"] = int(codon.split(":")[0])
    variant_dict["mutation_type"] = mutation_type
    variant_dict["name"] = name
    variant_dict["codon"] = ""
    variant_dict["mutation"] = mutation
    variant_dict["length"] = insdel_length
    variant_dict["rejected"] = rejected

    if rejected == 1:
        logging.info(
            "Rejected insdel with del len %s, ins len %s, and insdel len %s at pos %s with variant %s.",
            deletion_length,
            insertion_length,
            insdel_length,
            pos,
            variant,
        )

    return variant_dict


def process_variants_file(gatk_list, designed_variants_df, ref_AA_sequence, max_deletion_length):
    #   -Reject any variants with multiple substitutions
    #   -Reject any frameshifting mutations
    #   (designed mutations are singles, interpreting multi-codon indels as singles)
    #   -Accept arbitrary in-frame insertions
    #   -Accept arbitrary in-frame deletions
    #   -Further process "insdel" calls
    # For each one, we will add the observed counts
    # to the appropriate
    # For synonymous variants, the "length_codon" is 0

    # input:
    #          gatk_list: the gatk file, as a list of lists (one line per element)
    #
    #           designed_variants_df: the expected variants dataframe
    #               fields:
    #                 int count - number of observations of variant
    #                   int pos - position of variant (in codons)
    #                   str mutation_type - mutation class: S = synonymous
    #                                                       M = substitution
    #                                                       I = insertion
    #                                                       D = deletion
    #                 str name - mutation name, in 1 letter amino acids
    #                 str codon - for substitutions, the designed codon. empty for others
    #                 str mutation - for substitutions, the variant residue. for
    #                                indels, the specific variant, in the form (I/D)_(1/2/3),
    #                                as appropriate
    #                 int length - length of change, in codons
    #                 str hgvs - hgvs string of mutation
    #
    #           ref_sequence: the reference amino acid sequence, as a string
    #
    #           max_deletion_length: the maximum length of a deletion to be considered

    # output: variants_df with fields as above, with counts added to each variant
    #         rejected_list - list of rejected variants. Contents of each line are
    #                         the gatk output lines.
    #         counts_stats - dict of counts statistics

    # Keep track of counts processed and rejected for statistics

    rejected_stats = {}
    accepted_stats = {}
    total_stats = {}

    rejected_stats["outside_orf_counts"] = 0
    rejected_stats["fs_counts"] = 0
    rejected_stats["wrong_codon_counts"] = 0
    rejected_stats["wrong_variant_counts"] = 0
    rejected_stats["insdel_variant_counts"] = 0
    rejected_stats["multi_variant_counts"] = 0

    accepted_stats["accepted_syn_counts"] = 0
    accepted_stats["accepted_sub_counts"] = 0
    accepted_stats["accepted_stop_counts"] = 0
    accepted_stats["accepted_ins_counts"] = 0
    accepted_stats["accepted_del_counts"] = 0
    accepted_stats["accepted_insdel_counts"] = 0

    total_stats["total_counts"] = 0
    total_stats["total_rejected_counts"] = 0
    total_stats["total_accepted_counts"] = 0

    rejected_list = []
    variants_df = designed_variants_df.copy(deep=True)

    variant_names_array = sorted(variants_df["name"].array)

    for line in gatk_list:
        counts = int(line[0])
        length_NT = int(line[3])
        NT = line[4]
        length_codon = int(line[5])
        codon = line[6]
        AA = line[7]
        mutation = line[8]

        length = -1
        mutation_type = ""
        variant = ""
        name = ""
        count = 0
        pos = -1
        rejected = 1

        total_stats["total_counts"] = total_stats["total_counts"] + counts

        # If a variant is mapped outside of the ORF, nucleotide change length is 0 and
        # codon/AA fields are empty

        if length_NT == 0:
            rejected_list = rejected_list + [line]
            rejected_stats["outside_orf_counts"] = (
                rejected_stats["outside_orf_counts"] + counts
            )
            continue
        if AA == "":
            rejected_list = rejected_list + [line]
            rejected_stats["outside_orf_counts"] = (
                rejected_stats["outside_orf_counts"] + counts
            )
            continue

        # If a variant includes a frame shift, do not keep it.

        if "FS" in AA:
            rejected_stats["fs_counts"] = rejected_stats["fs_counts"] + counts
            rejected_list = rejected_list + [line]
            continue

        # Consider how many discontiguous non-synonymous mutations are found. If more
        # than 1, then this is not a designed variant. If none are found, this is a
        # synonymous variant.

        if mutation == "":
            # synonymous check
            if AA[0] == "S":
                count = counts
                pos = int(codon.split(":")[0])
                length = 1
                mutation_type = "S"
                variant = AA[-1]
                name = variant + str(pos) + variant
                rejected = 0
                accepted_stats["accepted_syn_counts"] = (
                    accepted_stats["accepted_syn_counts"] + counts
                )

        if len(mutation.split(";")) == 1:
            # Nonsyn check. Classify mutations here.
            # Note: checks for single position variants.

            # Is the variant a deletion?
            if mutation[-3:] == "del":
                deletion_dict = process_deletion(line)
                count = deletion_dict["counts"]
                pos = deletion_dict["pos"]
                length = deletion_dict["length"]
                mutation_type = deletion_dict["mutation_type"]
                variant = deletion_dict["mutation"]
                name = deletion_dict["name"]
                rejected = 0
                accepted_stats["accepted_del_counts"] = (
                    accepted_stats["accepted_del_counts"] + counts
                )

            # Is the variant an insertion?

            ins_re = regex.search(
                r"([a-zA-Z])([0-9]+)_([a-zA-Z])([0-9]+)ins([A-Z]+)", mutation
            )

            if ins_re:
                insertion_dict = process_insertion(line)
                count = insertion_dict["counts"]
                pos = insertion_dict["pos"]
                length = insertion_dict["length"]
                mutation_type = insertion_dict["mutation_type"]
                variant = insertion_dict["mutation"]
                name = insertion_dict["name"]
                rejected = 0
                accepted_stats["accepted_ins_counts"] = (
                    accepted_stats["accepted_ins_counts"] + counts
                )

            # Is the variant a substitution?

            if AA[0] == "M":
                if len(AA.split(", ")) == 1:
                    count = counts
                    pos = int(codon.split(":")[0])
                    length = 1
                    mutation_type = "M"
                    variant = AA[-1]
                    name = mutation
                    rejected = 0
                    accepted_stats["accepted_sub_counts"] = (
                        accepted_stats["accepted_sub_counts"] + counts
                    )

            # Is the variant a stop?

            if AA[0] == "N":
                if len(AA.split(", ")) == 1:
                    count = counts
                    pos = int(codon.split(":")[0])
                    length = 1
                    mutation_type = "N"
                    variant = AA[-1]
                    name = mutation
                    rejected = 0
                    accepted_stats["accepted_stop_counts"] = (
                        accepted_stats["accepted_stop_counts"] + counts
                    )

            # # Check if the mutation is a designed one before adding.

            if mutation_type == "M" or mutation_type == "S":
                try:
                    designed_codon = variants_df[
                        (variants_df["pos"] == pos)
                        & (variants_df["mutation"] == variant)
                    ]["codon"].array[0]
                    # If this codon wasn't one we designed, don't add it (that is, add 0)
                    if designed_codon != codon[-3:]:
                        rejected_stats["wrong_codon_counts"] = (
                            rejected_stats["wrong_codon_counts"] + counts
                        )
                        rejected = 1
                        count = 0
                        # If this codon was one we designed, but the mutation type is wrong, remove counts from stats
                        if mutation_type == "M":
                            accepted_stats["accepted_sub_counts"] = (
                                accepted_stats["accepted_sub_counts"] - counts
                            )
                        else:
                            accepted_stats["accepted_syn_counts"] = (
                                accepted_stats["accepted_syn_counts"] - counts
                            )
                except IndexError:
                    rejected_stats["wrong_variant_counts"] = (
                        rejected_stats["wrong_variant_counts"] + counts
                    )
                    count = 0
                    rejected = 1
                    if mutation_type == "M":
                        accepted_stats["accepted_sub_counts"] = (
                            accepted_stats["accepted_sub_counts"] - counts
                        )
                    else:
                        accepted_stats["accepted_syn_counts"] = (
                            accepted_stats["accepted_syn_counts"] - counts
                        )

            insdel_re = regex.search(
                r"([a-zA-Z])([0-9]+)_([a-zA-Z])([0-9]+)insdel([a-zA-Z]+)", mutation
            )

            if insdel_re:
                insdel_dict = process_insdel(line, variant_names_array, ref_AA_sequence, max_deletion_length)
                count = insdel_dict["counts"]
                pos = insdel_dict["pos"]
                length = insdel_dict["length"]
                mutation_type = insdel_dict["mutation_type"]
                variant = insdel_dict["mutation"]
                name = insdel_dict["name"]
                rejected = insdel_dict["rejected"]

                if rejected == 1:
                    rejected_stats["insdel_variant_counts"] = (
                        rejected_stats["insdel_variant_counts"] + counts
                    )
                else:
                    accepted_stats["accepted_insdel_counts"] = (
                        accepted_stats["accepted_insdel_counts"] + counts
                    )

            # Add counts to the variant df, if they are not flagged for rejection
            # No need to add counts here, since they should already by added for rejected ones.

            if rejected:
                rejected_list = rejected_list + [line]
                continue

            try:
                variants_df.loc[variants_df.name == name, "count"] += int(count)
            except:
                if not rejected:
                    rejected = 1
                    rejected_stats["wrong_variant_counts"] = rejected_stats[
                        "wrong_variant_counts"
                    ] + int(count)
                else:
                    pass

        # If a multiple-change variant
        else:
            rejected_list = rejected_list + [line]
            rejected_stats["multi_variant_counts"] = (
                rejected_stats["multi_variant_counts"] + counts
            )

    total_stats["total_rejected_counts"] = (
        rejected_stats["outside_orf_counts"]
        + rejected_stats["fs_counts"]
        + rejected_stats["wrong_codon_counts"]
        + rejected_stats["wrong_variant_counts"]
        + rejected_stats["insdel_variant_counts"]
        + rejected_stats["multi_variant_counts"]
    )
    total_stats["total_accepted_counts"] = (
        accepted_stats["accepted_syn_counts"]
        + accepted_stats["accepted_sub_counts"]
        + accepted_stats["accepted_stop_counts"]
        + accepted_stats["accepted_ins_counts"]
        + accepted_stats["accepted_del_counts"]
        + accepted_stats["accepted_insdel_counts"]
    )

    return variants_df, rejected_list, rejected_stats, accepted_stats, total_stats


def write_enrich_df(file, variant_df):
    wt_summed = variant_df[variant_df["mutation_type"] == "S"]["count"].sum()

    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        variant_df.to_csv(f, columns=["hgvs", "count"], index=False, sep="\t")
        f.write("_wt\t" + str(wt_summed))

    return


def write_stats_file(file, stats_dict):
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        for key, value in stats_dict.items():
            f.write(key + "\t" + str(value) + "\n")

    return
