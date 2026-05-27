import csv
import logging
import os
import tempfile
import unittest
from unittest.mock import patch

# Import the functions you want to test directly
# Adjust if the script is named differently or in a different location
from workflow.rules.scripts.generate_variants import (
    designed_variants,
    extract_codon,
    get_sequence_segment,
    name_to_hgvs,
    write_designed_csv,
)


class TestOligoProcessing(unittest.TestCase):
    def setUp(self):
        # Set up common test data
        self.reference_sequence = (
            "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        )
        self.offset = 1  # 1-based offset for testing

    def test_name_to_hgvs(self):
        """Test conversion of variant name to HGVS format."""
        self.assertEqual(name_to_hgvs("A1B"), "p.(A1B)")
        self.assertEqual(name_to_hgvs("M154L"), "p.(M154L)")
        self.assertEqual(name_to_hgvs("K23_Y24del"), "p.(K23_Y24del)")

    def test_get_sequence_segment_linear(self):
        """Test getting sequence segments for linear genomes."""
        ref = "ABCDEFGHIJ"

        # Normal case
        self.assertEqual(get_sequence_segment(ref, 2, 5, False), "CDE")

        # Edge cases
        self.assertEqual(get_sequence_segment(ref, 0, 3, False), "ABC")
        self.assertEqual(get_sequence_segment(ref, 7, 10, False), "HIJ")

        # Out of bounds
        self.assertEqual(get_sequence_segment(ref, -2, 3, False), "ABC")
        self.assertEqual(get_sequence_segment(ref, 8, 15, False), "IJ")

    def test_get_sequence_segment_circular(self):
        """Test getting sequence segments for circular genomes."""
        ref = "ABCDEFGHIJ"

        # Normal case (no wrap)
        self.assertEqual(get_sequence_segment(ref, 2, 5, True), "CDE")

        # Wrap around the origin
        self.assertEqual(get_sequence_segment(ref, 8, 12, True), "IJAB")
        self.assertEqual(get_sequence_segment(ref, -2, 3, True), "IJABC")

        # Full wrap around
        self.assertEqual(get_sequence_segment(ref, 8, 18, True), "IJABCDEFGH")

    def test_extract_codon(self):
        """Test extracting codons from oligo sequences."""
        # Test case where pre_codon is used
        oligo = "GACTCGATCGATATGCTAGCA"
        pre_codon = "GACTCGATCGA"
        self.assertEqual(extract_codon(oligo, pre_codon, "", 0, "", 0, 1), "TAT")

        # Test case where post_codon is used
        oligo = "ATGCTAGCAGACTCGATCGA"
        post_codon = "GACTCGATCGA"
        self.assertEqual(extract_codon(oligo, "", post_codon, 0, "", 0, 1), "GCA")

        # If both flanking sequences are empty, the function returns early.
        ref = "ATGCTAGCAGACTCGATCGA"
        oligo = "XXATGCTAGCAXXXX"
        self.assertEqual(extract_codon(oligo, "", "", 0, ref, 0, 1, False), "")

        # Test with no matches
        self.assertEqual(extract_codon("ACTG", "XXX", "YYY", 0, "ZZZZ", 0, 1), "")

    def test_extract_codon_flanks_with_regex_metacharacters(self):
        """Flanking sequences must be escaped before re.split.

        A pre_codon containing a regex metacharacter (e.g. `*`, `.`, `(`)
        would either crash with re.error or match the wrong bases. The
        characters here aren't legal IUPAC, but a user-supplied reference
        FASTA can plausibly contain any of them (alignment gaps, stop-codon
        annotations, ambiguity notations like `(A/G)`), so the function
        must not assume clean ACGT input.
        """
        oligo = "AAA.BBBTATCCC"
        # `.` is the regex wildcard; without escaping, `re.split("AAA.BBB", ...)`
        # would split anywhere AAA + any-char + BBB appears, which here happens
        # to match but in general gives the wrong answer or fails to anchor.
        pre_codon = "AAA.BBB"
        self.assertEqual(
            extract_codon(oligo, pre_codon, "", 0, "", 0, 1),
            "TAT",
        )

        # `(` is unbalanced and would raise re.error("missing ), unterminated subpattern")
        # without escaping.
        oligo_with_paren = "AAA(GROUP)TATCCC"
        self.assertEqual(
            extract_codon(oligo_with_paren, "AAA(GROUP)", "", 0, "", 0, 1),
            "TAT",
        )

    def test_designed_variants_substitution(self):
        """Test processing substitution variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                b"name,sequence\ntest_DMS-1_Ala10Gly,ACTAGCTAGCGCTAGCTAGCT\n"
            )

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        # Call the function
        with patch(
            "workflow.rules.scripts.generate_variants.extract_codon",
            return_value="GGC",
        ):
            variants = designed_variants(temp_name, ref, offset)

        # Verify the results
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0]["mutation_type"], "M")  # Missense
        self.assertEqual(variants[0]["name"], "A10G")
        self.assertEqual(variants[0]["codon"], "GGC")
        self.assertEqual(variants[0]["hgvs"], "p.(A10G)")

        os.unlink(temp_name)

    def test_designed_variants_deletion(self):
        """Test processing deletion variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                b"name,sequence\ntest_delete-1_3-5,ACTAGCTAGCGCTAGCTAGCT\n"
            )

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        variants = designed_variants(temp_name, ref, offset)

        # Verify the results
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0]["mutation_type"], "D")  # Deletion
        self.assertEqual(variants[0]["pos"], 5)
        self.assertEqual(variants[0]["mutant"], "D_1")  # 1 codon deletion
        self.assertEqual(variants[0]["length"], 1)

        os.unlink(temp_name)

    def test_designed_variants_deletion_missing_position_is_skipped(self):
        """Malformed deletion oligo (missing trailing -POSITION) must not crash
        the pipeline. DIMPLE always emits all three components, so we treat the
        two-component form as a non-match rather than raising on int(None)."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            # Note: no trailing "-<pos>" — would previously crash with
            # TypeError: int() argument must be a string ... not 'NoneType'
            temp_file.write(
                b"name,sequence\ntest_delete-1_3,ACTAGCTAGCGCTAGCTAGCT\n"
            )

        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        # Should return cleanly with no matched variant, not raise.
        variants = designed_variants(temp_name, ref, offset)
        self.assertEqual(variants, [])

        os.unlink(temp_name)

    def test_designed_variants_insertion(self):
        """Test processing insertion variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                b"name,sequence\ntest_insert-1_ATGGCG-5,ACTAGCTAGCGCTAGCTAGCT\n"
            )

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        variants = designed_variants(temp_name, ref, offset)

        # Verify the results
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0]["mutation_type"], "I")  # Insertion
        self.assertTrue("ins" in variants[0]["name"])
        self.assertEqual(variants[0]["codon"], "ATGGCG")
        self.assertEqual(variants[0]["mutant"], "I_2")  # 2 codon insertion

        os.unlink(temp_name)

    # ------------------------------------------------------------------
    # Insertion length-consistency check.
    #
    # The oligo name declares the inserted sequence; the actual oligo
    # nucleotide sequence is the ground truth. Anchored on flanking
    # reference codons, the bases between the flanks in the oligo are
    # the actual inserted region. A length mismatch indicates the
    # library's name and content disagree (replaces the prior dead
    # `if insertion_length != len(inserted_seq)` self-comparison).
    # ------------------------------------------------------------------

    # A reference with distinct codons so the script's flanking-window
    # anchors aren't ambiguous (the more-realistic-looking
    # "ATGGCTAGCATGGCTAGC..." reference used by the other tests has a
    # tandem repeat that makes pre_window match in multiple positions —
    # fine for those tests, but it would mask the length-check logic).
    _NONREPEATING_REF = "AAACCCGGGTTTAATCATCAGTACAAGTGGAGCAGTGAACGTTCC"

    def _build_insertion_oligo(self, ref, offset, pos, actual_inserted):
        """Construct an oligo whose actual inserted region equals
        `actual_inserted`, with flanking taken from `ref` around `pos`."""
        start_codon_pos = offset + (3 * (pos - 1))
        end_codon_pos = offset + (3 * pos)
        # PRE_SPAN is 15 in the script; we use the full available flank
        # so the script's re.split anchors find a unique match.
        pre = ref[max(0, start_codon_pos - 15) : start_codon_pos + 3]
        post = ref[end_codon_pos : end_codon_pos + 3 + 15]
        return pre + actual_inserted + post

    def test_designed_variants_insertion_length_mismatch_warns(self):
        """Oligo NAME declares a different number of inserted bases than the
        oligo SEQUENCE actually contains — emit a length-mismatch warning."""
        ref = self._NONREPEATING_REF
        offset = 1
        pos = 5

        # Oligo contains a 3-base insertion ("TGC") but the name declares
        # a 9-base inserted sequence ("AAATTTGGG").
        oligo_seq = self._build_insertion_oligo(ref, offset, pos, "TGC")
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                f"name,sequence\ntest_insert-1_AAATTTGGG-{pos},{oligo_seq}\n".encode()
            )

        with self.assertLogs(level="WARNING") as cm:
            designed_variants(temp_name, ref, offset)

        self.assertTrue(
            any("length mismatch" in msg.lower() for msg in cm.output),
            f"Expected an insertion length mismatch warning. Got: {cm.output}",
        )

        os.unlink(temp_name)

    def test_designed_variants_insertion_length_match_no_warning(self):
        """When the name and the oligo sequence agree, no mismatch warning."""
        ref = self._NONREPEATING_REF
        offset = 1
        pos = 5
        inserted = "TGCATG"  # 6 bases

        oligo_seq = self._build_insertion_oligo(ref, offset, pos, inserted)
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                f"name,sequence\ntest_insert-1_{inserted}-{pos},{oligo_seq}\n".encode()
            )

        with self.assertLogs(level="WARNING") as cm:
            designed_variants(temp_name, ref, offset)
            # Sentinel so assertLogs doesn't fail when the function is silent.
            logging.warning("__sentinel__")

        mismatch_warnings = [
            msg for msg in cm.output if "length mismatch" in msg.lower()
        ]
        self.assertEqual(
            mismatch_warnings,
            [],
            f"Unexpected length-mismatch warning: {mismatch_warnings}",
        )

        os.unlink(temp_name)

    def test_designed_variants_insertion_unfindable_flank_warns(self):
        """If the oligo doesn't contain the expected flanking ref codons,
        the script can't locate the insertion and must say so explicitly
        rather than silently skipping the length check."""
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1
        # Oligo sequence completely unrelated to the reference.
        oligo_seq = "NNNNNNNNNNNNNNNNNN"
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name
            temp_file.write(
                f"name,sequence\ntest_insert-1_ATGGCG-5,{oligo_seq}\n".encode()
            )

        with self.assertLogs(level="WARNING") as cm:
            designed_variants(temp_name, ref, offset)

        self.assertTrue(
            any("anchor" in msg.lower() and "flank" in msg.lower() for msg in cm.output),
            f"Expected a flank-anchoring warning. Got: {cm.output}",
        )

        os.unlink(temp_name)

    def test_write_designed_csv(self):
        """Test writing variants to a CSV file."""
        # Create sample variants
        variants = [
            [0, 10, "M", "A10G", "GGC", "G", 1, "p.(A10G)", 1],
            [0, 20, "D", "M20del", "", "D_1", 1, "p.(M20del)", 2],
        ]

        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name

        # Call the function
        header = [
            "count",
            "pos",
            "mutation_type",
            "name",
            "codon",
            "mutation",
            "length",
            "hgvs",
            "chunk",
        ]
        write_designed_csv(temp_name, header, variants)

        # Verify the output
        with open(temp_name, "r") as f:
            reader = csv.reader(f)
            rows = list(reader)

            # Check header
            self.assertEqual(rows[0], header)

            # Check data rows
            self.assertEqual(rows[1][0], "0")  # index
            self.assertEqual(rows[1][1], "10")  # pos
            self.assertEqual(rows[1][2], "M")  # mutation_type
            self.assertEqual(rows[1][3], "A10G")  # name

            self.assertEqual(rows[2][0], "0")  # index
            self.assertEqual(rows[2][1], "20")  # pos
            self.assertEqual(rows[2][2], "D")  # mutation_type
            self.assertEqual(rows[2][3], "M20del")  # name

        os.unlink(temp_name)


class TestIntegration(unittest.TestCase):
    """Integration tests to test the full workflow."""

    def setUp(self):
        # Create sample reference sequence
        self.ref_sequence = "ATGGCTAGCAAGTTAGCGCTGAGCTAG"

        # Create temporary files
        with tempfile.NamedTemporaryFile(delete=False) as ref_file:
            self.ref_file_name = ref_file.name
            ref_file.write(b">reference\n")
            ref_file.write(self.ref_sequence.encode())

        with tempfile.NamedTemporaryFile(delete=False) as oligo_file:
            self.oligo_file_name = oligo_file.name
            oligo_file.write(b"name,sequence\n")
            oligo_file.write(b"test_DMS-1_Ala3Gly,ATGGCTGGCAAGTTAGCGCTGAGCTAG\n")
            oligo_file.write(b"test_delete-1_3-5,ATGGCTAGCTAGCGCTGAGCTAG\n")
            oligo_file.write(
                b"test_insert-1_ATGGCG-7,ATGGCTAGCAAGTTAGCATGGCGCTGAGCTAG\n"
            )

        with tempfile.NamedTemporaryFile(delete=False) as output_file:
            self.output_file_name = output_file.name

    def tearDown(self):
        # Clean up temporary files
        os.unlink(self.ref_file_name)
        os.unlink(self.oligo_file_name)
        os.unlink(self.output_file_name)

    def test_full_workflow(self):
        """Test the full workflow with mocked snakemake object."""

        # This would normally execute the script, but we'll test the components separately
        # For a full integration test, you'd need to refactor the script to make it importable

        # For demonstration, we'll just run the main components manually
        variants = designed_variants(
            self.oligo_file_name, self.ref_sequence, 1, False  # offset  # is_circular
        )

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

        write_designed_csv(self.output_file_name, header, variants)

        # Verify output file exists
        self.assertTrue(os.path.exists(self.output_file_name))


if __name__ == "__main__":
    unittest.main()
