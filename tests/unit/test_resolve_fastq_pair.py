"""
Tests for resolve_fastq_pair — the helper that maps a sample-prefix string
to a paired-end R1/R2 fastq pair on disk.

The historical bug: an unanchored `data_dir.glob(f"{filename}*")` would
treat a prefix like `1_S1` as a partial match against `1_S10_*` /
`1_S100_*`, etc. Depending on which side picked up the wrong neighbor,
this either silently cross-paired (wrong R1 + wrong R2 for the same
nominal sample) or tripped the "Multiple R1 files found" error. The fix
anchors the glob on a `.` or `_` boundary immediately after the prefix.
"""

from pathlib import Path

import pytest

from workflow.rules.scripts.script_utils import resolve_fastq_pair


def _touch(dir: Path, *names: str) -> None:
    for n in names:
        (dir / n).touch()


def test_resolves_basic_illumina_pair(tmp_path):
    _touch(tmp_path, "sampleA_R1_001.fastq.gz", "sampleA_R2_001.fastq.gz")
    r1, r2 = resolve_fastq_pair(tmp_path, "sampleA")
    assert r1.name == "sampleA_R1_001.fastq.gz"
    assert r2.name == "sampleA_R2_001.fastq.gz"


def test_rejects_partial_prefix_match(tmp_path):
    """The reported regression: `1_S1` must not match `1_S10_*`."""
    _touch(
        tmp_path,
        "1_S10_R1_001.fastq.gz",
        "1_S10_R2_001.fastq.gz",
        "1_S100_R1_001.fastq.gz",
        "1_S100_R2_001.fastq.gz",
    )
    with pytest.raises(FileNotFoundError, match="prefix '1_S1'"):
        resolve_fastq_pair(tmp_path, "1_S1")


def test_resolves_correct_pair_when_longer_neighbors_exist(tmp_path):
    """Both `1_S1` and `1_S10` should resolve independently to their own
    files when both are present."""
    _touch(
        tmp_path,
        "1_S1_R1_001.fastq.gz",
        "1_S1_R2_001.fastq.gz",
        "1_S10_R1_001.fastq.gz",
        "1_S10_R2_001.fastq.gz",
    )
    r1, r2 = resolve_fastq_pair(tmp_path, "1_S1")
    assert r1.name == "1_S1_R1_001.fastq.gz"
    assert r2.name == "1_S1_R2_001.fastq.gz"

    r1, r2 = resolve_fastq_pair(tmp_path, "1_S10")
    assert r1.name == "1_S10_R1_001.fastq.gz"
    assert r2.name == "1_S10_R2_001.fastq.gz"


def test_handles_simplified_naming(tmp_path):
    """{prefix}_R1.fastq.gz (no _001 suffix)."""
    _touch(tmp_path, "expt_R1.fastq.gz", "expt_R2.fastq.gz")
    r1, r2 = resolve_fastq_pair(tmp_path, "expt")
    assert r1.name == "expt_R1.fastq.gz"
    assert r2.name == "expt_R2.fastq.gz"


def test_handles_numeric_naming(tmp_path):
    """{prefix}_1.fastq.gz / {prefix}_2.fastq.gz."""
    _touch(tmp_path, "expt_1.fastq.gz", "expt_2.fastq.gz")
    r1, r2 = resolve_fastq_pair(tmp_path, "expt")
    assert r1.name == "expt_1.fastq.gz"
    assert r2.name == "expt_2.fastq.gz"


def test_handles_dot_separator(tmp_path):
    """{prefix}.R1.fastq.gz instead of {prefix}_R1.fastq.gz."""
    _touch(tmp_path, "expt.R1.fastq.gz", "expt.R2.fastq.gz")
    r1, r2 = resolve_fastq_pair(tmp_path, "expt")
    assert r1.name == "expt.R1.fastq.gz"
    assert r2.name == "expt.R2.fastq.gz"


def test_handles_fq_gz_extension(tmp_path):
    _touch(tmp_path, "expt_R1.fq.gz", "expt_R2.fq.gz")
    r1, r2 = resolve_fastq_pair(tmp_path, "expt")
    assert r1.name == "expt_R1.fq.gz"
    assert r2.name == "expt_R2.fq.gz"


def test_handles_uncompressed_fastq(tmp_path):
    _touch(tmp_path, "expt_R1.fastq", "expt_R2.fastq")
    r1, r2 = resolve_fastq_pair(tmp_path, "expt")
    assert r1.name == "expt_R1.fastq"
    assert r2.name == "expt_R2.fastq"


def test_raises_when_no_files_match(tmp_path):
    with pytest.raises(FileNotFoundError, match="prefix 'missing'"):
        resolve_fastq_pair(tmp_path, "missing")


def test_raises_when_only_r1_present(tmp_path):
    _touch(tmp_path, "expt_R1.fastq.gz")
    with pytest.raises(FileNotFoundError):
        resolve_fastq_pair(tmp_path, "expt")


def test_raises_on_duplicate_r1(tmp_path):
    """If two files both look like R1 for the same prefix, that's an error."""
    _touch(
        tmp_path,
        "expt_R1.fastq.gz",
        "expt_R1_001.fastq.gz",
        "expt_R2.fastq.gz",
    )
    with pytest.raises(ValueError, match="Multiple R1"):
        resolve_fastq_pair(tmp_path, "expt")
