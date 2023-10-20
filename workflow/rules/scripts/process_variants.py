import csv
import pathlib

import pandas as pd
import regex

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


def process_variants_file(gatk_list, designed_variants_df):
    #   -Reject any variants with multiple substitutions
    #   -Reject any frameshifting mutations
    #   (designed mutations are singles, interpreting multi-codon indels as singles)
    #   -Accept G, GS, and GSG insertions
    #   -Accept 1,2,3x deletions
    #   -Further process "insdel" calls
    # For each one, we will add the observed counts
    # to the appropriate
    # For synonymous variants, the "length_codon" is 0

    # input: variants_df fields:
    #           int count - number of observations of variant
    #           int pos - position of variant (in codons)
    #           int chunk_pos - position of variant modulo "chunk"
    #           int chunk - chunk number of variant
    #           str mutation_type - mutation class: S = synonymous
    #                                      M = substitution
    #                                      I = insertion
    #                                      D = deletion
    #           str name - mutation name, in 1 letter amino acids
    #           str codon - for substitutions, the designed codon. empty for others
    #           str mutation - for substitutions, the variant residue. for
    #                           indels, the specific variant, in the form (I/D)_(1/2/3),
    #                           as appropriate
    #           int length - length of change, in codons
    #           str hgvs - hgvs string of mutation

    rejected_list = []
    variants_df = designed_variants_df.copy(deep=True)

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

        # If a variant is mapped outside of the ORF, nucleotide change length is 0 and
        # codon/AA fields are empty

        if length_NT == 0:
            pass

        # Consider how many discontiguous non-synonymous mutations are found. If more
        # than 1, then this is nota designed variant. If none are found, this is a
        # synonymous variant.

        if AA == "":
            continue

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

        if len(mutation.split(";")) == 1:
            # Nonsyn check. Classify mutations here.

            # Is the variant a deletion?
            if mutation[-3:] == "del":
                if length_codon <= 3:
                    count = counts
                    pos = int(codon.split(":")[0])
                    length = length_codon
                    mutation_type = "D"
                    variant = "D_" + str(length)
                    name = mutation
                    rejected = 0

            # Is the variant an insertion?
            # TODO: this doesn't need to check for AA[0] being "I", first.
            # But this should work for now.
            AAs = AA.split(",")
            for m in AAs:
                if m.strip()[0] == "I":
                    # Check for designed variants
                    if mutation[-4:] == "insG":
                        length = 1
                    if mutation[-5:] == "insGS":
                        length = 2
                    if mutation[-6:] == "insGSG":
                        length = 3

                    if length > 0:
                        rejected = 0
                        count = counts
                        pos = int(codon.split(":")[0])
                        mutation_type = "I"
                        variant = "I_" + str(length)
                        name = mutation

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

            # # Check if the nonsynonymous mutation is a designed one before adding.

            if mutation_type == "M" or mutation_type == "S":
                try:
                    designed_codon = variants_df[
                        (variants_df["pos"] == pos)
                        & (variants_df["mutation"] == variant)
                    ]["codon"].array[0]
                    # If this codon wasn't one we designed, don't add it (that is, add 0)
                    if designed_codon != codon[-3:]:
                        count = 0
                        rejected = 1
                except IndexError:
                    count = 0
                    rejected = 1

            # Is variant mapping having a hard time with insdel?
            # This accounts for the case where deletion of the
            # form XYZA->X--- (just deletion) is being called as
            # XYZA->---X (deletion + mutation)
            # Note: bbmap defaults to left alignment, as is standard for NGS reads.
            # This means we only need to consider incorrectly left-aligned cases.

            # TODO: consider case with insdelZZ or insdelZZZ

            insdel_re = regex.search(
                r"([a-zA-Z])([0-9]+)_([a-zA-Z])([0-9]+)insdel([a-zA-Z]+)", mutation
            )

            if insdel_re:
                count = 0

                length = int(insdel_re.group(4)) - int(insdel_re.group(2))
                insdel_length = len(insdel_re.group(5))

                if (
                    insdel_re.group(1) == insdel_re.group(5)
                    or insdel_re.group(3) == insdel_re.group(5)
                    or insdel_re.group(1) == insdel_re.group(3)
                ):
                    if length - insdel_length < 5:
                        rejected = 0
                        count = counts
                        # pos = int(codon.split(":")[0]) + 1
                        pos = int(insdel_re.group(2)) + insdel_length
                        length = length
                        mutation_type = "D"
                        variant = "D_" + str(length)
                        name = mutation

                        # Find the correct canonical name

                        try:
                            name = variants_df[
                                (variants_df["pos"] == pos)
                                & (variants_df["mutation"] == variant)
                            ]["name"].array[0]

                        except:
                            rejected = 1
                            # print("Name error:", name)

            # Add counts to the variant df

            try:
                variants_df.loc[variants_df.name == name, "count"] += int(count)
            except:
                rejected = 1

            if rejected:
                rejected_list = rejected_list + [line]

        # If a multiple-change variant
        else:
            rejected_list = rejected_list + [line]

    return variants_df, rejected_list


def write_enrich_df(file, variant_df):
    wt_summed = variant_df[variant_df["mutation_type"] == "S"]["count"].sum()

    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        variant_df.to_csv(f, columns=["hgvs", "count"], index=False, sep="\t")
        f.write("_wt\t" + str(wt_summed))

    return


def remove_zeros_enrich(enrich_file_list):
    # Input is a list of files containing Enrich2 formatted tsvs.
    # This function combines them all into a single df, then removes the variants
    # that are missing across all experiments, and rewrites the tsvs.
    # Warning: this function operates IN-PLACE!

    df_list = [
        pd.read_csv(f, sep="\t", names=["hgvs", f], header=0).set_index("hgvs") for f in enrich_file_list
    ]

    combined_enrich_df = pd.concat(df_list, join="outer", axis=1)

    unobserved_variants = combined_enrich_df[
        (combined_enrich_df == 0).all(axis=1)
    ].index.to_list()

    combined_enrich_df = combined_enrich_df[(combined_enrich_df != 0).any(axis=1)]

    # write out each file again

    for enrich_file in enrich_file_list:
        p = pathlib.Path(enrich_file)
        p.parent.mkdir(parents=True, exist_ok=True)

        with p.open("w+") as f:
            combined_enrich_df.to_csv(
                f, columns=[enrich_file], header=["count"], index=True, sep="\t"
            )

    return unobserved_variants
