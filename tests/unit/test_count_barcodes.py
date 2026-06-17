"""Unit tests for the barcode-counting core (count_barcodes.py).

Covers the pure pieces of barcode mode: synonymous-HGVS recognition, the two
extractor flavors (regex / fixed-position), FASTQ streaming (plain + gzip),
read tallying with the unknown/no-tag/min-count dispositions, variant
aggregation over the full map universe, and the mutation_type synthesis that
feeds write_enrich_df's _wt row.
"""

import gzip

import pandas as pd
import process_variants

from workflow.rules.scripts.count_barcodes import (
    aggregate_variant_counts,
    build_counts_frame,
    count_barcodes,
    is_synonymous,
    iter_fastq_sequences,
    make_extractor,
)


def _write_fastq(path, seqs):
    """Write a minimal FASTQ (gzip if .gz) with the given sequences."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


class TestIsSynonymous:
    def test_population_level(self):
        assert is_synonymous("p.(=)")

    def test_explicit_equals(self):
        assert is_synonymous("p.(L54=)")

    def test_same_aa_substitution(self):
        assert is_synonymous("p.(L54L)")

    def test_missense_is_not(self):
        assert not is_synonymous("p.(L54P)")

    def test_deletion_is_not(self):
        assert not is_synonymous("p.(M1del)")


class TestMakeExtractor:
    def test_regex_captures_group(self):
        extract = make_extractor(barcode_pattern="^(.{4})")
        assert extract("AAAACCCCGG") == "AAAA"

    def test_regex_miss_returns_none(self):
        extract = make_extractor(barcode_pattern="^(ACGT)")
        assert extract("TTTTTTTT") is None

    def test_fixed_position_slice(self):
        extract = make_extractor(barcode_start=2, barcode_length=3)
        assert extract("AAAACCCC") == "AAC"

    def test_fixed_position_too_short_returns_none(self):
        extract = make_extractor(barcode_start=4, barcode_length=4)
        assert extract("AAA") is None


class TestIterFastqSequences:
    def test_plain_fastq(self, tmp_path):
        path = _write_fastq(tmp_path / "reads.fastq", ["AAAA", "CCCC"])
        assert list(iter_fastq_sequences(path)) == ["AAAA", "CCCC"]

    def test_gzipped_fastq_uppercased(self, tmp_path):
        path = _write_fastq(tmp_path / "reads.fastq.gz", ["aaaa", "cccc"])
        assert list(iter_fastq_sequences(path)) == ["AAAA", "CCCC"]


class TestCountBarcodes:
    MAP = {"AAAA": "p.(G10A)", "CCCC": "p.(R11K)"}

    def _extract(self):
        return make_extractor(barcode_start=0, barcode_length=4)

    def test_matching_and_dispositions(self, tmp_path):
        path = _write_fastq(
            tmp_path / "r.fastq.gz",
            ["AAAATTTT", "AAAATTTT", "CCCCGGGG", "GGGGAAAA", "AA"],
        )
        counts, stats = count_barcodes(path, self.MAP, self._extract())
        assert counts == {"AAAA": 2, "CCCC": 1}
        assert stats["total_reads"] == 5
        assert stats["matched_reads"] == 3
        assert stats["unknown_barcode"] == 1  # GGGG not in map
        assert stats["no_barcode_extracted"] == 1  # "AA" too short
        assert stats["distinct_barcodes_observed"] == 2

    def test_min_barcode_count_floor(self, tmp_path):
        path = _write_fastq(
            tmp_path / "r.fastq.gz", ["AAAA", "AAAA", "AAAA", "CCCC"]
        )
        counts, stats = count_barcodes(
            path, self.MAP, self._extract(), min_barcode_count=2
        )
        assert counts == {"AAAA": 3}  # CCCC seen once, below floor -> dropped
        assert stats["below_min_count_barcodes"] == 1
        assert stats["below_min_count_reads"] == 1


class TestAggregateVariantCounts:
    def test_full_universe_with_zeros(self):
        # Two barcodes map to the same variant; a third variant is unobserved.
        mapping = {"AAAA": "p.(G10A)", "TTTT": "p.(G10A)", "CCCC": "p.(R11K)"}
        barcode_counts = {"AAAA": 5, "TTTT": 3}
        variant_counts = aggregate_variant_counts(barcode_counts, mapping)
        assert variant_counts == {"p.(G10A)": 8, "p.(R11K)": 0}


class TestBuildCountsFrame:
    def test_mutation_type_synthesized(self):
        df = build_counts_frame({"p.(L54L)": 4, "p.(L54P)": 2})
        row = df.set_index("hgvs")
        assert row.loc["p.(L54L)", "mutation_type"] == "S"
        assert row.loc["p.(L54P)", "mutation_type"] == "M"
        assert set(df.columns) == {"hgvs", "count", "mutation_type"}
        assert isinstance(df, pd.DataFrame)


class TestEnrichWriterIntegration:
    """build_counts_frame -> process_variants.write_enrich_df is the seam where
    barcode counts join the existing scoring path; verify the emitted TSV has
    the (hgvs, count) rows plus the synonymous-sum _wt row scoring expects."""

    def test_enrich_tsv_has_counts_and_wt_row(self, tmp_path):
        df = build_counts_frame({"p.(G10A)": 7, "p.(L54L)": 3, "p.(L54=)": 2})
        out = tmp_path / "sample.tsv"
        process_variants.write_enrich_df(out, df, noprocess=False)

        text = out.read_text()
        lines = text.strip().splitlines()
        assert lines[0] == "hgvs\tcount"
        body = dict(line.split("\t") for line in lines[1:])
        assert body["p.(G10A)"] == "7"
        # _wt is the sum of the two synonymous variants (3 + 2).
        assert body["_wt"] == "5"
