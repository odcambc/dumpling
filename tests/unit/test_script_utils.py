import pytest

from workflow.rules.scripts.script_utils import file_digest, run_script


def test_run_script_logs_and_reraises(mock_snakemake, mocker):
    snakemake = mock_snakemake(log=["/dev/null"])
    basic_config = mocker.patch("workflow.rules.scripts.script_utils.logging.basicConfig")
    log_exception = mocker.patch("workflow.rules.scripts.script_utils.logging.exception")

    def boom(_snakemake):
        raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        run_script(snakemake, boom)

    basic_config.assert_called_once()
    log_exception.assert_called_once_with("Script failed")


class TestFileDigest:
    """Hash helper used to key the bbmap index path on reference content."""

    def test_same_content_same_digest(self, tmp_path):
        a = tmp_path / "a.fasta"
        b = tmp_path / "b.fasta"
        content = b">ref\nACGTACGTACGT\n"
        a.write_bytes(content)
        b.write_bytes(content)
        assert file_digest(a) == file_digest(b)

    def test_different_content_different_digest(self, tmp_path):
        a = tmp_path / "ref_v1.fasta"
        b = tmp_path / "ref_v2.fasta"
        a.write_bytes(b">ref\nACGTACGTACGT\n")
        b.write_bytes(b">ref\nACGTACGTACGA\n")  # one base different
        assert file_digest(a) != file_digest(b)

    def test_default_length_is_12_hex_chars(self, tmp_path):
        f = tmp_path / "ref.fasta"
        f.write_bytes(b">ref\nACGT\n")
        digest = file_digest(f)
        assert len(digest) == 12
        assert all(c in "0123456789abcdef" for c in digest), digest

    def test_length_kwarg_overrides_default(self, tmp_path):
        f = tmp_path / "ref.fasta"
        f.write_bytes(b">ref\nACGT\n")
        assert len(file_digest(f, length=8)) == 8
        assert len(file_digest(f, length=40)) == 40  # full sha1

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            file_digest(tmp_path / "does_not_exist.fasta")

    def test_streams_in_chunks_for_large_files(self, tmp_path):
        # Should not load the whole file into memory; correctness check
        # is the same digest as if computed in one shot.
        import hashlib

        f = tmp_path / "big.fasta"
        # 1 MB of varied bytes — well above the 64 KB read chunk in file_digest.
        content = (b"ACGTN" * 209715)[:1_048_576]
        f.write_bytes(content)
        expected = hashlib.sha1(content).hexdigest()[:12]
        assert file_digest(f) == expected
