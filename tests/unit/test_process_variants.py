# test_process_variants.py

import pytest
import pandas as pd

from workflow.rules.scripts.process_variants import (
    read_gatk_csv,
    process_variants_file,
    parse_multi_synonymous,
    process_insertion,
    process_deletion,
    process_insdel,
    process_single_site,
    write_enrich_df,
    write_stats_file,
    name_to_hgvs,
)

# -------------------------
#     FIXTURES
# -------------------------


@pytest.fixture
def designed_variants_df():
    """
    A small mock DataFrame representing the "designed variants" in your pipeline.
    Each row corresponds to a single variant the pipeline expects.
    """
    data = {
        "count": [0, 0, 0],
        "pos": [10, 11, 12],
        "mutation_type": ["M", "S", "D"],
        "name": ["G10A", "R11R", "P12del"],
        "codon": ["AAA", "CGT", ""],
        "mutation": ["A", "R", "D_1"],
        "length": [1, 1, 1],
        "hgvs": ["p.(G10A)", "p.(R11R)", "p.(P12del)"],
    }
    df = pd.DataFrame(data)
    return df


@pytest.fixture
def ref_aa_sequence():
    """
    Example reference AA sequence. In real usage,
    you'd read from a FASTA or config.
    """
    return "MKVAFWLLLS"


# -------------------------
#     TESTS
# -------------------------


def test_name_to_hgvs():
    """Check the name_to_hgvs function's string formatting."""
    assert name_to_hgvs("A12B") == "p.(A12B)"
    assert name_to_hgvs("") == "p.()"
    assert name_to_hgvs("X999X") == "p.(X999X)"


def test_parse_multi_synonymous():
    """Check that parse_multi_synonymous properly builds the multi-syn name from codons."""
    codon_string = "97:CAG>CAC, 98:GGA>GGG"
    result = parse_multi_synonymous(codon_string)
    # Each codon is translated to single-letter AA
    # For example:
    #  CAG -> Q
    #  CAC -> H
    #  GGA -> G
    #  GGG -> G
    # So we'd get "Q97H;G98G" or something similar
    # We'll just check the structure.
    # The exact translation depends on your code, but here's an example expectation:
    assert ";" in result, "Expected multiple codons separated by semicolons"
    # We can do a rough check:
    parts = result.split(";")
    assert len(parts) == 2, "Should have two parsed codons"
    assert "97" in parts[0]
    assert "98" in parts[1]


def test_process_insertion():
    """Check that an insertion is parsed to the correct dictionary structure."""
    # A mock GATK line. We'll only set the relevant columns
    # Based on your variantCounts_colnames:
    line = [
        "10",  # counts
        "0",  # coverage
        "0",  # mean_length
        "0",  # length_NT
        "A",  # NT
        "2",  # length_codon
        "10:AAA>CCC, 11:TTT>GGG",  # codon
        "Mstuff",  # AA
        "M10_11insCCC",  # mutations
    ]
    result = process_insertion(line)
    assert result["mutation_type"] == "I"
    assert result["count"] == 10
    assert result["pos"] == 10
    assert result["mutation"] == "I_2"


def test_process_deletion():
    """Check that a deletion is parsed to the correct dictionary structure."""
    line = [
        "5",  # counts
        "0",
        "0",
        "0",
        "A",
        "1",  # length_codon
        "10:AAA>---",  # codon
        "Mstuff",
        "A10del",
    ]
    result = process_deletion(line)
    assert result["mutation_type"] == "D"
    assert result["count"] == 5
    assert result["length"] == 1
    assert result["rejected"] is False


def test_process_insdel():
    """Check that an insdel is handled. Note that your logic is specialized."""
    line = [
        "7",
        "0",
        "0",
        "0",
        "A",
        "2",  # length_codon
        "10:AAA>TTT",  # codon
        "FSstuff",  # AA (maybe a frameshift?)
        "A10_B11insdelCC",
    ]
    # For demonstration, we pass noprocess=True
    result = process_insdel(line, "MAGICREFSEQ", noprocess=True)
    assert result["count"] == 7
    assert result["mutation_type"] in ["ID", "D", "Z"], "Depending on your logic"
    assert "insdel" in result["mutation"] or "ID_" in result["mutation"]
    assert result["rejected"] is False


def test_process_single_site_no_mutation(ref_aa_sequence):
    """Check single site logic for a no-mutation (synonymous) scenario."""
    line = [
        "3",  # counts
        "0",
        "0",
        "3",
        "A",
        "1",
        "10:AAA>AAA",
        "Sstuff",  # means 'synonymous' if it starts with S
        "",  # mutation is empty
    ]
    out = process_single_site(line, ref_aa_sequence, noprocess=False)
    assert out["mutation_type"] == "S"
    assert out["hgvs"].startswith("p.(")


def test_process_variants_file_noprocess(designed_variants_df, ref_aa_sequence):
    """Check that process_variants_file, with noprocess=True, accumulates counts
    into a new DataFrame without filtering based on designed_variants_df."""
    # We'll create a mock GATK list with multiple lines
    mock_gatk_list = [
        ["5", "0", "0", "3", "A", "1", "10:AAA>AAA", "M:G>A", "G10A"],
        ["2", "0", "0", "3", "A", "1", "11:AAA>AAA", "S:R>R", ""],  # synonymous
        [
            "10",
            "0",
            "0",
            "3",
            "A",
            "1",
            "13:AAA>---",
            "D:F>-",
            "F13del",
        ],
    ]
    # This last line might be 'rejected' under normal process, but noprocess should keep it
    variants_df, rejected_list, rejected_stats, accepted_stats, total_stats = (
        process_variants_file(
            mock_gatk_list,
            designed_variants_df,
            ref_aa_sequence,
            max_deletion_length=3,
            noprocess=True,
        )
    )

    # We expect the final variants_df to have 3 rows
    assert len(variants_df) == 3
    # Check total counts
    assert total_stats["total_counts"] == 17  # 5 + 2 + 10
    # Under noprocess, we normally don't do rejections, except for orf/FS checks
    # The 3rd line has FSstuff in 'AA', but your code might or might not skip it
    # because we do minimal checks. Adjust as needed.


def test_process_variants_file_normal(designed_variants_df, ref_aa_sequence):
    """Check that process_variants_file, with noprocess=False, accumulates counts only
    for 'expected' variants in designed_variants_df."""
    mock_gatk_list = [
        ["5", "0", "0.1", "1", "1774:A>C", "1", "10:CAA>AAG", "M:G>A", "G10A"],
        [
            "2",
            "0",
            "1.2",
            "2",
            "1774:A>C, 1823:G>A",
            "1",
            "11:CGC>CGT",
            "S:R>R",
            "",
        ],  # synonymous
        [
            "10",
            "0",
            "0",
            "3",
            "1774:A>-, 1775:A>-, 1776:A>-",
            "1",
            "13:AAA>---",
            "D:F>-",
            "F13del",
        ],
    ]
    variants_df, rejected_list, rejected_stats, accepted_stats, total_stats = (
        process_variants_file(
            mock_gatk_list,
            designed_variants_df,
            ref_aa_sequence,
            max_deletion_length=3,
            noprocess=False,
        )
    )

    # We expect the first line to add 5 to the "G10A" row
    row_G10A = variants_df.loc[variants_df["name"] == "G10A"]
    assert row_G10A["count"].iloc[0] == 5

    # Third line is not in designed_variants_df, so it should be rejected
    assert len(rejected_list) == 1

    # Check stats
    assert total_stats["total_counts"] == 17  # 5 + 2 + 10
    # You can check sub-counters in accepted_stats / rejected_stats if needed:
    # e.g., accepted_stats["accepted_sub_counts"] == 5, etc.


# -------------------------
#  H2: vectorized accumulator regression coverage
# -------------------------


def _gatk_g10a(counts):
    """Build a GATK row that should be recognized as the G10A variant."""
    return [str(counts), "0", "0.1", "1", "1774:A>C", "1", "10:CAA>AAG", "M:G>A", "G10A"]


def test_process_variants_file_sums_duplicate_observations(
    designed_variants_df, ref_aa_sequence
):
    """The hot path used to do `variants_df.loc[mask, "count"] += counts` once
    per GATK row. The vectorized rewrite accumulates into a dict and merges
    once at the end. This guards the invariant that K duplicate observations
    of the same variant produce a count of sum(K), not the most recent K."""
    gatk_list = [_gatk_g10a(3), _gatk_g10a(7), _gatk_g10a(2)]
    variants_df, *_ = process_variants_file(
        gatk_list,
        designed_variants_df,
        ref_aa_sequence,
        max_deletion_length=3,
        noprocess=False,
    )
    g10a = variants_df.loc[variants_df["name"] == "G10A", "count"].iloc[0]
    assert g10a == 12  # 3 + 7 + 2


def test_process_variants_file_unobserved_variants_have_count_zero(
    designed_variants_df, ref_aa_sequence
):
    """A designed variant that's never observed must keep its initial count
    of 0 — the vectorized merge should not perturb rows whose names don't
    appear in the GATK list."""
    gatk_list = [_gatk_g10a(5)]
    variants_df, *_ = process_variants_file(
        gatk_list,
        designed_variants_df,
        ref_aa_sequence,
        max_deletion_length=3,
        noprocess=False,
    )
    # G10A observed, R11R / P12del not observed.
    assert variants_df.loc[variants_df["name"] == "G10A", "count"].iloc[0] == 5
    assert variants_df.loc[variants_df["name"] == "R11R", "count"].iloc[0] == 0
    assert variants_df.loc[variants_df["name"] == "P12del", "count"].iloc[0] == 0


def test_process_single_site_unexpected_aa_is_dropped(ref_aa_sequence):
    """A row with mutation='' and AA[0] not in {S, M, N} previously fell through
    process_single_site with count=0, silently discarding the observation. M10:
    such rows are now flagged with mutation_type='X', rejected=True, and the
    original count is preserved so the main loop can route them to a dedicated
    rejected-stats bucket."""
    line = [
        "11",  # counts
        "0", "0", "3", "A", "1",
        "10:AAA>TTT",
        "Xstuff",  # AA[0] is 'X', not in {S, M, N}
        "",  # mutation empty
    ]
    out = process_single_site(line, ref_aa_sequence, noprocess=False)
    assert out["mutation_type"] == "X"
    assert out["rejected"] is True
    assert out["count"] == 11


def test_process_variants_file_routes_unexpected_to_dedicated_stat(
    designed_variants_df, ref_aa_sequence
):
    """The fall-through case (mutation='' and AA[0] not in {S,M,N}) must be
    accounted for in rejected_stats['unexpected_mutation_counts'] and rolled
    into total_rejected_counts, not silently absorbed into wrong_variant_counts
    or — worse — dropped from the totals entirely. Realistic input: a normal
    synonymous row alongside two unexpected rows."""
    mock_gatk_list = [
        ["3", "0", "0", "3", "A", "1", "11:AAA>AAA", "S:R>R", ""],
        ["4", "0", "0", "3", "A", "1", "10:AAA>TTT", "Xstuff", ""],
        ["7", "0", "0", "3", "A", "1", "11:AAA>CCC", "?bogus", ""],
    ]
    variants_df, rejected_list, rejected_stats, _, total_stats = process_variants_file(
        mock_gatk_list,
        designed_variants_df,
        ref_aa_sequence,
        max_deletion_length=3,
        noprocess=True,
    )
    assert rejected_stats["unexpected_mutation_counts"] == 11  # 4 + 7
    assert len(rejected_list) == 2
    assert total_stats["total_counts"] == 14  # 3 + 4 + 7
    assert total_stats["total_rejected_counts"] == 11
    # The dropped rows must not appear in the noprocess output dataframe;
    # the synonymous row survives.
    assert len(variants_df) == 1


def test_process_variants_file_all_rejected_under_noprocess_does_not_crash(
    designed_variants_df, ref_aa_sequence
):
    """Defensive: when every GATK row is rejected under noprocess=True the
    output dataframe is empty. The post-loop column cleanup must tolerate
    that — previously it raised KeyError because `.drop(columns=['rejected'])`
    has nothing to drop on a no-columns frame."""
    mock_gatk_list = [
        ["4", "0", "0", "3", "A", "1", "10:AAA>TTT", "Xstuff", ""],
    ]
    variants_df, *_ = process_variants_file(
        mock_gatk_list,
        designed_variants_df,
        ref_aa_sequence,
        max_deletion_length=3,
        noprocess=True,
    )
    assert len(variants_df) == 0


def test_process_variants_file_order_independent(
    designed_variants_df, ref_aa_sequence
):
    """Shuffling the GATK row order must not change the final counts."""
    rows = [_gatk_g10a(c) for c in (1, 5, 9, 2, 4)]

    def run(gatk_list):
        df, *_ = process_variants_file(
            gatk_list,
            designed_variants_df,
            ref_aa_sequence,
            max_deletion_length=3,
            noprocess=False,
        )
        return df.loc[df["name"] == "G10A", "count"].iloc[0]

    forward = run(rows)
    reverse = run(list(reversed(rows)))
    shuffled = run([rows[2], rows[0], rows[4], rows[1], rows[3]])
    assert forward == reverse == shuffled == 21  # 1 + 5 + 9 + 2 + 4
