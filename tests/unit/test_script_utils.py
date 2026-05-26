import pandas as pd
import pytest
from Bio.Seq import Seq

from workflow.rules.scripts.script_utils import (
    file_digest,
    load_experiments,
    run_script,
    translate_orf,
    validate_experiment_time_or_bin,
    validate_scoring_backend_mode,
)


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


class TestTranslateOrf:
    """One translation helper, shared between generate_variants and
    process_counts, so the two scripts can't disagree on what the ORF's
    AA sequence is — especially for circular ORFs that cross the origin."""

    def test_linear_orf_translates_inclusive_1_based_slice(self):
        # GCT|TTT|... in DNA -> A|F|... in protein
        ref = Seq("GCTTTT")
        assert str(translate_orf(ref, "1-6")) == "AF"

    def test_linear_orf_subrange(self):
        # Skip the first codon: TTT -> F only.
        ref = Seq("GCTTTT")
        assert str(translate_orf(ref, "4-6")) == "F"

    def test_circular_orf_wraps_origin(self):
        # 6 nt ring; ORF from position 5 to position 4 wraps: indexes
        # 4..5 (CC) + 0..3 (GCTT) = "CCGCTT" = P|L when translated.
        ref = Seq("GCTTCC")
        # generate_variants's circular concatenation:
        # ref[5-1:] + ref[:4] = "CC" + "GCTT" = "CCGCTT" -> "PL"
        assert str(translate_orf(ref, "5-4")) == "PL"

    def test_circular_and_linear_diverge_when_wrapping(self):
        # The whole point: a linear slice with start > end is empty in
        # Python, but the circular interpretation produces the real ORF.
        # The pre-fix process_counts.py used to return "" for this case
        # while generate_variants returned the right answer; both should
        # now agree.
        ref = Seq("GCTTCC")
        linear_naive = ref[5 - 1 : 4].translate()
        circular = translate_orf(ref, "5-4")
        assert str(linear_naive) == ""
        assert str(circular) != ""


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


class TestLoadExperiments:
    """Canonical experiment-CSV loader shared by every rule that touches the file."""

    def test_indexed_by_sample_with_sample_column_retained(self, tmp_path):
        csv = tmp_path / "exp.csv"
        csv.write_text("sample,condition,time,file\ns1,c,0,f1\ns2,c,1,f2\n")
        df = load_experiments(csv)
        assert list(df.index) == ["s1", "s2"]
        assert "sample" in df.columns

    def test_strips_utf8_bom(self, tmp_path):
        """Excel-exported CSVs ship with a BOM; the loader must transparently
        eat it so downstream columns are not renamed to '\\ufeffsample'."""
        csv = tmp_path / "exp_bom.csv"
        csv.write_bytes(
            b"\xef\xbb\xbfsample,condition,time,file\ns1,c,0,f1\n"
        )
        df = load_experiments(csv)
        assert "sample" in df.columns
        assert df.loc["s1", "condition"] == "c"

    def test_drops_blank_rows(self, tmp_path):
        csv = tmp_path / "exp_blanks.csv"
        csv.write_text(
            "sample,condition,time,file\ns1,c,0,f1\n,,,\ns2,c,1,f2\n,,,\n"
        )
        df = load_experiments(csv)
        assert list(df.index) == ["s1", "s2"]

    def test_duplicate_sample_raises(self, tmp_path):
        csv = tmp_path / "exp_dup.csv"
        csv.write_text("sample,condition,time,file\ns1,c,0,f1\ns1,c,1,f2\n")
        with pytest.raises(ValueError, match="Duplicate values found in 'sample'"):
            load_experiments(csv)


class TestValidateExperimentTimeOrBin:
    """Replaces the schema-level `oneOf` so error messages name the offending sample."""

    def test_pure_timecourse_passes(self):
        df = pd.DataFrame({
            "sample": ["a", "b"],
            "time": [0, 1],
        })
        validate_experiment_time_or_bin(df)

    def test_pure_facs_passes(self):
        df = pd.DataFrame({
            "sample": ["a", "b"],
            "bin": [1, 2],
        })
        validate_experiment_time_or_bin(df)

    def test_mixed_rows_pass_when_each_row_picks_one(self):
        df = pd.DataFrame({
            "sample": ["timecourse_row", "facs_row"],
            "time": [0, None],
            "bin": [None, 2],
        })
        validate_experiment_time_or_bin(df)

    def test_row_with_neither_fails_with_sample_name(self):
        df = pd.DataFrame({
            "sample": ["s_ok", "s_bad"],
            "time": [0, None],
        })
        with pytest.raises(ValueError, match="'s_bad'.*missing 'time' or 'bin'"):
            validate_experiment_time_or_bin(df)

    def test_row_with_both_fails_with_sample_name(self):
        df = pd.DataFrame({
            "sample": ["s_ok", "s_both"],
            "time": [0, 5],
            "bin": [None, 3],
        })
        with pytest.raises(ValueError, match="'s_both'.*both set"):
            validate_experiment_time_or_bin(df)

    def test_missing_columns_entirely_fails_listing_every_sample(self):
        df = pd.DataFrame({
            "sample": ["s1", "s2", "s3"],
            "condition": ["c", "c", "c"],
        })
        with pytest.raises(ValueError) as excinfo:
            validate_experiment_time_or_bin(df)
        msg = str(excinfo.value)
        for s in ["s1", "s2", "s3"]:
            assert f"'{s}'" in msg

    def test_time_zero_is_valid_not_treated_as_missing(self):
        """time=0 is the standard T0 timepoint, not a null value."""
        df = pd.DataFrame({
            "sample": ["a"],
            "time": [0],
        })
        validate_experiment_time_or_bin(df)

    def test_all_problems_reported_at_once(self):
        """User shouldn't have to fix one row, re-run, fix the next."""
        df = pd.DataFrame({
            "sample": ["bad1", "ok", "bad2"],
            "time": [None, 0, 5],
            "bin": [None, None, 3],
        })
        with pytest.raises(ValueError) as excinfo:
            validate_experiment_time_or_bin(df)
        msg = str(excinfo.value)
        assert "'bad1'" in msg
        assert "'bad2'" in msg
        assert "'ok'" not in msg


class TestValidateScoringBackendMode:
    """Lilace's Stan model requires a synonymous control set and has no
    total-counts fallback. noprocess mode doesn't produce trustworthy
    synonymous labels, so the combination must be rejected at parse time.

    rosace-aa hits the same problem plus a stricter version: its
    RunRosace.Rosace signature requires wt.col/mut.col/ctrl.col +
    aa.code. The noprocess path produces none of those, so the combination
    is structurally invalid.
    """

    def test_lilace_with_noprocess_rejected(self):
        config = {"scoring_backend": "lilace", "noprocess": True}
        with pytest.raises(ValueError, match="lilace.*incompatible.*noprocess"):
            validate_scoring_backend_mode(config)

    def test_lilace_without_noprocess_passes(self):
        config = {"scoring_backend": "lilace", "noprocess": False}
        validate_scoring_backend_mode(config)

    def test_rosace_aa_with_noprocess_rejected(self):
        config = {"scoring_backend": "rosace_aa", "noprocess": True}
        with pytest.raises(ValueError, match="rosace_aa.*incompatible.*noprocess"):
            validate_scoring_backend_mode(config)

    def test_rosace_aa_without_noprocess_passes(self):
        config = {"scoring_backend": "rosace_aa", "noprocess": False}
        validate_scoring_backend_mode(config)

    def test_rosace_with_noprocess_passes(self):
        config = {"scoring_backend": "rosace", "noprocess": True}
        validate_scoring_backend_mode(config)

    def test_rosace_without_noprocess_passes(self):
        config = {"scoring_backend": "rosace", "noprocess": False}
        validate_scoring_backend_mode(config)

    def test_missing_keys_treated_as_falsy(self):
        """Defaults aren't this function's job — schema validation fills
        them in. An empty config shouldn't crash here."""
        validate_scoring_backend_mode({})
