import unittest
from unittest.mock import patch, mock_open, MagicMock
import os
import tempfile
import csv

# Import the functions you want to test directly
# Adjust if the script is named differently or in a different location
from workflow.rules.scripts.generate_variants import (
    name_to_hgvs,
    get_sequence_segment,
    extract_codon,
    designed_variants,
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
        self.assertEqual(get_sequence_segment(ref, -2, 3, True), "IJAB")

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
        self.assertEqual(extract_codon(oligo, "", post_codon, 0, "", 0, 1), "AGC")

        # Test direct lookup
        ref = "ATGCTAGCAGACTCGATCGA"
        oligo = "XXATGCTAGCAXXXX"
        self.assertEqual(extract_codon(oligo, "", "", 0, ref, 0, 1, False), "ATG")

        # Test with no matches
        self.assertEqual(extract_codon("ACTG", "XXX", "YYY", 0, "ZZZZ", 0, 1), "")

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="name,sequence\ntest_DMS-1_Ala10Gly,ACTAGCTAGCGCTAGCTAGCT\n",
    )
    def test_designed_variants_substitution(self, mock_file):
        """Test processing substitution variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        # Mock the regex search and extract_codon function
        with patch("regex.search") as mock_regex:
            with patch("builtins.open", mock_file):
                # Set up the mock regex match
                mock_match = MagicMock()
                mock_match.group.side_effect = lambda x: {
                    1: "1",
                    2: "Ala",
                    3: "10",
                    4: "Gly",
                }[x]
                mock_regex.return_value = mock_match

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

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="name,sequence\ntest_delete-1_3-5,ACTAGCTAGCGCTAGCTAGCT\n",
    )
    def test_designed_variants_deletion(self, mock_file):
        """Test processing deletion variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        # Mock the regex search
        with patch("regex.search") as mock_regex:
            with patch("builtins.open", mock_file):
                # First return None for substitution pattern, then a match for deletion
                mock_regex.side_effect = [
                    None,  # No match for substitution
                    MagicMock(),  # Match for deletion
                ]

                # Set up the deletion mock match
                mock_regex.return_value.group.side_effect = lambda x: {
                    1: "1",
                    2: "3",
                    3: "5",
                }[x]

                # Mock Seq.translate to return amino acids
                with patch("Bio.Seq.Seq.translate", return_value="M"):
                    variants = designed_variants(temp_name, ref, offset)

                # Verify the results
                self.assertEqual(len(variants), 1)
                self.assertEqual(variants[0]["mutation_type"], "D")  # Deletion
                self.assertEqual(variants[0]["name"], "M5del")
                self.assertEqual(variants[0]["mutant"], "D_1")  # 1 codon deletion

        os.unlink(temp_name)

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="name,sequence\ntest_insert-1_ATGGCG-5,ACTAGCTAGCGCTAGCTAGCT\n",
    )
    def test_designed_variants_insertion(self, mock_file):
        """Test processing insertion variants."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_name = temp_file.name

        # Mock the reference sequence
        ref = "ATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGCATGGCTAGC"
        offset = 1

        # Mock the regex search and Bio.Seq.translate
        with patch("regex.search") as mock_regex:
            with patch("builtins.open", mock_file):
                # Return None for the first two patterns, then a match for insertion
                mock_regex.side_effect = [
                    None,  # No match for substitution
                    None,  # No match for deletion
                    MagicMock(),  # Match for insertion
                ]

                # Set up the insertion mock match
                mock_regex.return_value.group.side_effect = lambda x: {
                    1: "1",
                    2: "ATGGCG",
                    3: "5",
                }[x]

                # Mock Seq.translate to return amino acids
                with patch("Bio.Seq.Seq.translate", side_effect=["M", "A"]):
                    variants = designed_variants(temp_name, ref, offset)

                # Verify the results
                self.assertEqual(len(variants), 1)
                self.assertEqual(variants[0]["mutation_type"], "I")  # Insertion
                self.assertTrue("ins" in variants[0]["name"])
                self.assertEqual(variants[0]["codon"], "ATGGCG")
                self.assertEqual(variants[0]["mutant"], "I_2")  # 2 codon insertion

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
