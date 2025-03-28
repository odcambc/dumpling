from typing import Dict, List, Tuple, Any, Union
import csv
import pathlib
import logging

import pandas as pd
import regex

import Bio.Seq

# Type aliases
VariantDict = Dict[str, Union[int, str, bool]]
StatsDict = Dict[str, int]
GatkList = List[List[str]]

aa_3to1_dict: Dict[str, str] = {
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
aa_1to3_dict: Dict[str, str] = {v: k for k, v in aa_3to1_dict.items()}

variantCounts_colnames: List[str] = [
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


def name_to_hgvs(name: str) -> str:
    return "p.(" + name + ")"


def read_gatk_csv(file: Union[str, pathlib.Path]) -> GatkList:
    # GATK output unfortunately may not be full-width, which will confuse pandas during
    # reading in. We'll read each line and make sure it's full-length before passing it
    # to pandas.
    logging.info("Loading GATK output file: %s", file)

    gatk_list: GatkList = []

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


def process_insertion(line: List[str]) -> VariantDict:
    # This contains the logic for parsing insertion variants of arbitrary length.
    # We define the position of an insertion as the position of the first new codon
    # in the new frame.

    variant_dict: VariantDict = {}

    count = int(line[0])
    length = int(line[5])
    mutation_type = "I"
    variant = line[8]
    mutation = "I_" + str(length)
    codons = line[6]
    pos = int(codons.split(":")[0])
    rejected = False

    codon = ""

    if length > 1:
        for c in codons.split(", "):
            codon = codon + c.split(">")[1]
    else:
        codon = codon + codons.split(">")[1]

    variant_dict["count"] = count
    variant_dict["pos"] = pos
    variant_dict["mutation_type"] = mutation_type
    variant_dict["name"] = variant
    variant_dict["codon"] = codon
    variant_dict["mutation"] = mutation
    variant_dict["length"] = length
    variant_dict["rejected"] = rejected

    return variant_dict


def process_deletion(line: List[str]) -> VariantDict:
    # This contains the logic for parsing deletion variants of arbitrary length.
    # We define the position of deletion as the position of the first deleted codon.

    variant_dict: VariantDict = {}
    codon = line[6]
    length_codon = int(line[5])
    variant = line[8]
    rejected = False

    if "_" in variant:
        start_codon = int(variant.split("_")[0][1:])
        end_codon = int(variant.split("_")[1][1:-3])

        # It may be the case that deletions are coupled with synonymous variants.
        # Reject these.
        if length_codon != (end_codon - start_codon + 1):
            rejected = True
    elif length_codon > 1:
        rejected = True

    variant_dict["count"] = int(line[0])
    variant_dict["pos"] = int(codon.split(":")[0])
    variant_dict["mutation_type"] = "D"
    variant_dict["name"] = variant
    variant_dict["codon"] = ""
    variant_dict["mutation"] = "D_" + str(length_codon)
    variant_dict["length"] = length_codon
    variant_dict["rejected"] = rejected

    return variant_dict


def process_insdel(
    line: List[str], ref_AA_sequence: str, noprocess: bool
) -> VariantDict:
    # This contains the logic for parsing insdel variants and checking for an edge
    # case where in-frame deletions are being called as insdels.
    # Deletions of the form XYZA->X--- (just deletion) may be called as
    # XYZA->---X (deletion + mutation)
    # Note: bbmap defaults to left alignment, as is standard for NGS reads.
    # This means we only need to consider incorrectly left-aligned cases.

    # The logic here requires passing an array of the designed variant names to
    # search for the proper canonical name. It also requires returning an additional
    # "rejected" value in the dict.

    variant_dict: VariantDict = {}
    codon = line[6]
    length_codon = int(line[5])
    variant = line[8]
    count = int(line[0])
    rejected = True
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
    if not insdel_re:
        logging.warning("Error in insdel parsing")
        logging.warning(line)
        return variant_dict

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
        rejected = False
        pos = start_pos + insertion_length
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

    else:
        if noprocess:
            pos = start_pos
            name = variant
            mutation_type = "ID"
            mutation = "ID_" + str(insdel_length)

    rejected = False

    variant_dict["count"] = count
    variant_dict["pos"] = int(codon.split(":")[0])
    variant_dict["mutation_type"] = mutation_type
    variant_dict["name"] = name
    variant_dict["codon"] = ""
    variant_dict["mutation"] = mutation
    variant_dict["length"] = insdel_length
    variant_dict["rejected"] = rejected

    if name == "":
        logging.warning(f"Error in insdel name calculation: HGVS consequence {line[8]}")
        logging.warning(f"Issue with line: {line}")

    return variant_dict


def process_single_site(
    line: List[str], ref_AA_sequence: str, noprocess: bool
) -> VariantDict:
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
    rejected = False

    hgvs = name_to_hgvs(mutation)

    if mutation[-3:] == "del":
        deletion_dict = process_deletion(line)
        deletion_dict["hgvs"] = hgvs

        return deletion_dict

    ins_re = regex.search(r"([a-zA-Z])([0-9]+)_([a-zA-Z])([0-9]+)ins([A-Z]+)", mutation)

    if ins_re:
        insertion_dict = process_insertion(line)
        insertion_dict["hgvs"] = hgvs

        return insertion_dict

    if "insdel" in mutation:
        insdel_dict = process_insdel(line, ref_AA_sequence, noprocess)
        insdel_dict["hgvs"] = hgvs

        return insdel_dict

    if mutation == "":
        if AA[0] == "S":
            synonymous_dict = {
                "count": counts,
                "pos": int(codon.split(":")[0]),
                "mutation_type": "S",
                "name": "",
                "codon": codon,
                "mutation": mutation,
                "length": len(AA.split(", ")),
                "rejected": False,
                "hgvs": "",
            }
            if len(AA.split(", ")) > 1:
                synonymous_dict["name"] = parse_multi_synonymous(codon)
            else:
                variant = AA[-1]
                synonymous_dict["name"] = variant + codon.split(":")[0] + variant
            synonymous_dict["hgvs"] = name_to_hgvs(synonymous_dict["name"])

            return synonymous_dict

    codon_length = len(codon.split(", "))
    if codon_length > 1:
        AAs = AA.split(", ")
        for i in range(len(AAs)):
            if AAs[i][0] != "S":
                AA = AAs[i]
                codon = codon.split(", ")[i]
                pos = int(codon.split(":")[0])
                break

    if AA[0] == "M":
        if len(AA.split(", ")) == 1:
            count = counts
            pos = int(codon.split(":")[0])
            length = 1
            mutation_type = "M"
            variant = AA[-1]
            name = mutation
            rejected = False

    if AA[0] == "N":
        if len(AA.split(", ")) == 1:
            count = counts
            pos = int(codon.split(":")[0])
            length = 1
            mutation_type = "N"
            variant = AA[-1]
            name = mutation
            rejected = False

    if name == "":
        logging.warning("Error in variant parsing")
        logging.warning(line)

    return {
        "count": count,
        "pos": pos,
        "mutation_type": mutation_type,
        "name": name,
        "codon": codon,
        "mutation": mutation,
        "length": length,
        "rejected": rejected,
        "hgvs": hgvs,
    }


def parse_multi_synonymous(codon: str) -> str:
    codons = codon.split(", ")

    names = []

    for codon in codons:
        split_codon = codon.split(":")
        codon_position = int(split_codon[0])
        codon_aa_start = split_codon[1].split(">")[0]
        codon_aa_end = split_codon[1].split(">")[1]
        aa_start = Bio.Seq.translate(codon_aa_start)
        aa_end = Bio.Seq.translate(codon_aa_end)
        hgvs_string = aa_start + str(codon_position) + aa_end
        names = names + [hgvs_string]

    name = ";".join(names)

    return name


def check_expected(variant_dict: VariantDict, variant_names_array: List[str]) -> bool:
    if variant_dict["name"] in variant_names_array:
        return False
    else:
        return True


def update_stats(
    accepted_stats: StatsDict, rejected_stats: StatsDict, variant_dict: VariantDict
) -> Tuple[StatsDict, StatsDict]:
    if variant_dict["rejected"]:
        rejected_stats["wrong_variant_counts"] = rejected_stats[
            "wrong_variant_counts"
        ] + int(variant_dict["count"])

    if variant_dict["mutation_type"] == "S":
        accepted_stats["accepted_syn_counts"] = accepted_stats[
            "accepted_syn_counts"
        ] + int(variant_dict["count"])
    elif variant_dict["mutation_type"] == "M":
        accepted_stats["accepted_sub_counts"] = accepted_stats[
            "accepted_sub_counts"
        ] + int(variant_dict["count"])
    elif variant_dict["mutation_type"] == "N":
        accepted_stats["accepted_stop_counts"] = accepted_stats[
            "accepted_stop_counts"
        ] + int(variant_dict["count"])
    elif variant_dict["mutation_type"] == "I":
        accepted_stats["accepted_ins_counts"] = accepted_stats[
            "accepted_ins_counts"
        ] + int(variant_dict["count"])
    elif variant_dict["mutation_type"] == "D":
        accepted_stats["accepted_del_counts"] = accepted_stats[
            "accepted_del_counts"
        ] + int(variant_dict["count"])
    elif variant_dict["mutation_type"] == "Z":
        accepted_stats["accepted_insdel_counts"] = accepted_stats[
            "accepted_insdel_counts"
        ] + int(variant_dict["count"])

    return accepted_stats, rejected_stats


def process_variants_file(
    gatk_list: GatkList,
    designed_variants_df: pd.DataFrame,
    ref_AA_sequence: str,
    max_deletion_length: int,
    noprocess: bool = True,
) -> Tuple[pd.DataFrame, List[List[str]], StatsDict, StatsDict, StatsDict]:
    rejected_stats: StatsDict = {}
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

    accepted_stats: StatsDict = {}
    total_stats: StatsDict = {}

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

    rejected_list: List[List[str]] = []

    # If we are not processing the variants, we just add all the counts to a new one.
    # Otherwise we will populate a copy of the designed variants dataframe.
    variants_noprocess_dict: Dict[str, VariantDict] = {}

    if noprocess:
        variants_df = pd.DataFrame()
        variant_names_array = []
    else:
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
        rejected = False

        total_stats["total_counts"] = total_stats["total_counts"] + counts

        if length_NT == 0 or AA == "":
            rejected_list = rejected_list + [line]
            rejected_stats["outside_orf_counts"] = (
                rejected_stats["outside_orf_counts"] + counts
            )
            rejected = True
            continue

        if "FS" in AA:
            rejected_stats["fs_counts"] = rejected_stats["fs_counts"] + counts
            rejected_list = rejected_list + [line]
            rejected = True
            continue

        if len(mutation.split(";")) == 1:
            variant_dict = process_single_site(line, ref_AA_sequence, noprocess)
        else:
            variant_dict = {}
            variant_dict["count"] = counts
            variant_dict["pos"] = -1
            variant_dict["mutation_type"] = "Z"
            variant_dict["name"] = mutation
            variant_dict["codon"] = codon.replace(",", ";")
            variant_dict["mutation"] = mutation
            variant_dict["length"] = -1
            variant_dict["hgvs"] = name_to_hgvs(mutation)
            variant_dict["rejected"] = False

            if variant_dict["mutation"] == "":
                variant_dict["mutation"] = "S"
                variant_dict["name"] = parse_multi_synonymous(variant_dict["codon"])
                if variant_dict["name"] == "":
                    logging.warning("Error in multi-synonymous parsing")
                    logging.warning(line)

        # If we aren't filtering variants, just add them to the dict.
        if noprocess:
            try:
                if variant_dict["name"] in variants_noprocess_dict:
                    variants_noprocess_dict[variant_dict["name"]]["count"] = int(
                        variant_dict["count"]
                    ) + int(variants_noprocess_dict[variant_dict["name"]]["count"])
                else:
                    variants_noprocess_dict[str(variant_dict["name"])] = variant_dict
            except KeyError:
                logging.warning("Error in variant processing")
                logging.warning(line)

        # If we are filtering the variants, we need to check if the variant is in the
        # designed variants dataframe. If it is not, we will reject it.
        else:
            variant_dict["rejected"] = check_expected(variant_dict, variant_names_array)

            accepted_stats, rejected_stats = update_stats(
                accepted_stats, rejected_stats, variant_dict
            )

            if variant_dict["rejected"]:
                rejected_list = rejected_list + [line]
                rejected_stats["wrong_variant_counts"] = (
                    rejected_stats["wrong_variant_counts"] + counts
                )
                continue

            try:
                variants_df.loc[
                    variants_df.name == variant_dict["name"], "count"
                ] += int(counts)
            except KeyError:
                rejected = True
                rejected_list = rejected_list + [line]
                rejected_stats["wrong_variant_counts"] = (
                    rejected_stats["wrong_variant_counts"] + counts
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

    if noprocess:
        variants_df = pd.DataFrame.from_dict(variants_noprocess_dict, orient="index")
        # Drop the "rejected" column
        variants_df.drop(columns=["rejected"], inplace=True)
        # Counts are ints.
        variants_df["count"] = variants_df["count"].astype(int)
        # rename mutation column to mutant
        variants_df.rename(columns={"mutation": "mutant"}, inplace=True)

    return variants_df, rejected_list, rejected_stats, accepted_stats, total_stats


def write_enrich_df(
    file: Union[str, pathlib.Path], variant_df: pd.DataFrame, noprocess: bool = True
) -> None:
    if not noprocess:
        wt_summed = variant_df[variant_df["mutation_type"] == "S"]["count"].sum()
    else:
        wt_summed = None

    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        variant_df.to_csv(f, columns=["hgvs", "count"], index=False, sep="\t")
        if wt_summed:
            f.write("_wt\t" + str(wt_summed))

    return


def write_stats_file(file: Union[str, pathlib.Path], stats_dict: StatsDict) -> None:
    p = pathlib.Path(file)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w+") as f:
        for key, value in stats_dict.items():
            f.write(key + "\t" + str(value) + "\n")

    return
