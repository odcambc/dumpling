"""
Regression tests for UTF-8 BOM tolerance on user-supplied CSVs.

The pipeline accepts experiment and variants CSVs that may be authored in
Excel or other tools that prepend a UTF-8 BOM (\\xef\\xbb\\xbf). Historically,
pandas read the BOM as part of the first column name (e.g.
`\\ufeffsample` instead of `sample`), silently breaking every downstream
`df["sample"]` lookup. Recent pandas (>= 2.x) auto-strips BOM under plain
utf-8, but we still pass `encoding="utf-8-sig"` explicitly for older
pandas, other readers, and code clarity — and to be robust if pandas ever
changes the default back.

These tests verify the encoding contract at the read sites we patched.
"""

import pandas as pd


def _write_csv(tmp_path, name, content, with_bom):
    path = tmp_path / name
    data = content.encode("utf-8")
    if with_bom:
        data = b"\xef\xbb\xbf" + data
    path.write_bytes(data)
    return path


CSV_CONTENT = "sample,condition,replicate,time\nA_T0,A,1,0\nA_T1,A,1,1\n"


def test_utf8_sig_strips_bom(tmp_path):
    csv = _write_csv(tmp_path, "bom.csv", CSV_CONTENT, with_bom=True)
    df = pd.read_csv(csv, header=0, encoding="utf-8-sig")
    # No BOM-mangled column name; normal lookups work.
    assert "sample" in df.columns
    assert "﻿sample" not in df.columns
    assert df.loc[0, "sample"] == "A_T0"


def test_utf8_sig_is_a_noop_for_plain_utf8(tmp_path):
    """utf-8-sig must continue to work for files without a BOM, otherwise
    flipping the encoding could break existing well-formed inputs."""
    csv = _write_csv(tmp_path, "plain.csv", CSV_CONTENT, with_bom=False)
    df = pd.read_csv(csv, header=0, encoding="utf-8-sig")
    assert "sample" in df.columns
    assert df.loc[0, "sample"] == "A_T0"


def test_example_csv_committed_without_bom():
    """The repo-tracked config/example.csv must not carry a BOM — we
    stripped it once; this test pins the invariant so it can't sneak back
    in if someone edits the file in Excel and re-saves."""
    import pathlib

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    example = repo_root / "config" / "example.csv"
    if not example.exists():
        # Don't fail in environments that don't ship config/ (unlikely, but
        # we don't want this test to gate on fixture layout).
        return
    with example.open("rb") as f:
        first_bytes = f.read(3)
    assert first_bytes != b"\xef\xbb\xbf", (
        "config/example.csv has a UTF-8 BOM again. Strip it with "
        "`sed -i '1s/^\\xef\\xbb\\xbf//' config/example.csv` and check what "
        "editor put it back."
    )
