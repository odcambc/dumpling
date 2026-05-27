"""
Behavioral tests for the FASTA-header normalization in the `prepare_reference`
rule (workflow/rules/ref.smk). The rule rewrites the first `>` line of a
user-supplied FASTA to `>{reference_name}` so the contig name matches the
filename stem — required by BBMap index build + GATK SAM/ref contig matching
(Upstream #19).

These tests exercise the rule's awk command directly rather than spinning up
snakemake, so they run in milliseconds. If the rule's shell template ever
gets refactored, the awk string below must move in lockstep — there's no
clean way to import a shell template into Python. The test docstrings and a
comment in ref.smk reference each other so the drift is easy to catch.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

# The exact awk command embedded in workflow/rules/ref.smk's prepare_reference
# rule. If you edit one, edit the other; the tests in this module are the
# tripwire if they fall out of sync.
NORMALIZE_AWK = (
    "awk -v name='{name}' "
    "'NR==1 && /^>/ {{ sub(/^>.*/, \">\" name); }} {{ print }}' "
    "{input}"
)


def _run_normalize(input_path: Path, reference_name: str) -> str:
    """Invoke the exact awk command from the prepare_reference rule and
    return its stdout (the normalized FASTA content)."""
    if shutil.which("awk") is None:
        pytest.skip("awk not available on this system")
    cmd = NORMALIZE_AWK.format(name=reference_name, input=str(input_path))
    # Snakemake's shell directive uses bash; match that for double-brace escape
    # behavior (the rule's `{{ ... }}` becomes literal `{ ... }` in shell).
    cmd = cmd.replace("{{", "{").replace("}}", "}")
    result = subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True
    )
    return result.stdout


class TestHeaderRewrite:
    def test_mismatched_header_rewritten(self, tmp_path):
        """The headline Upstream #19 case: header `>some_weird_contig` and
        filename `my_ref.fasta` → normalized output's first line is
        `>my_ref`."""
        src = tmp_path / "my_ref.fasta"
        src.write_text(">some_weird_contig_name_42\nACGTACGT\n")
        out = _run_normalize(src, "my_ref")
        lines = out.splitlines()
        assert lines[0] == ">my_ref"
        # Sequence MUST be preserved byte-for-byte — only the header changes.
        assert lines[1] == "ACGTACGT"

    def test_matching_header_unchanged(self, tmp_path):
        """Happy case: header already matches filename stem. Normalization
        should be a content-preserving no-op (idempotent)."""
        src = tmp_path / "example_ref.fasta"
        src.write_text(">example_ref\nACGTACGT\n")
        out = _run_normalize(src, "example_ref")
        assert out == ">example_ref\nACGTACGT\n"

    def test_multi_contig_only_first_rewritten(self, tmp_path):
        """DMS references are single-contig in practice, but multi-contig
        FASTAs must NOT have their secondary contigs silently renamed —
        that would break any downstream code that relies on the secondary
        contig names. Only line 1 is touched."""
        src = tmp_path / "my_ref.fasta"
        src.write_text(
            ">weird_first\n"
            "ACGTACGT\n"
            ">secondary_contig\n"
            "TTTTAAAA\n"
        )
        out = _run_normalize(src, "my_ref")
        assert out == (
            ">my_ref\n"
            "ACGTACGT\n"
            ">secondary_contig\n"
            "TTTTAAAA\n"
        )

    def test_sequence_byte_preserved(self, tmp_path):
        """Concordance gate: the sequence bytes must be byte-identical
        before and after normalization. Header changes only."""
        src = tmp_path / "ref.fasta"
        original_seq = "ACGTACGTACGTNNACGT" * 50
        src.write_text(f">whatever\n{original_seq}\n")
        out = _run_normalize(src, "ref")
        new_seq = "\n".join(out.splitlines()[1:])
        assert new_seq == original_seq

    def test_reference_name_with_special_chars(self, tmp_path):
        """Reference filenames in DMS data often include underscores,
        hyphens, digits (e.g., `SLCO1B1_ref`). The awk substitution must
        emit those literally without interpretation."""
        src = tmp_path / "SLCO1B1_ref.fasta"
        src.write_text(">weird\nACGT\n")
        out = _run_normalize(src, "SLCO1B1_ref")
        assert out.splitlines()[0] == ">SLCO1B1_ref"
