import hashlib
import logging
import re
from pathlib import Path

import pandas as pd

# Match the trailing R1/R2 marker in a fastq filename. Anchored at the right
# end of the filename so e.g. `_1` in `1_S1_R1_001.fastq.gz` doesn't trip
# the R1 detector before `_R1` does. `[._]` covers both Illumina-style
# `_R1_001.fastq.gz` and dot-separated `prefix.R1.fastq.gz`.
_FASTQ_R1_PATTERN = re.compile(r"[._](?:R1|1)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")
_FASTQ_R2_PATTERN = re.compile(r"[._](?:R2|2)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")


def resolve_fastq_pair(data_dir, filename):
    """Resolve a file prefix to a paired-end R1/R2 fastq pair.

    Handles common naming conventions:
      - Illumina standard: {prefix}_R1_001.fastq.gz
      - Simplified: {prefix}_R1.fastq.gz
      - Numeric: {prefix}_1.fastq.gz
      - Any of the above with .fq.gz, .fastq, or .fq extensions

    The glob is anchored on a `.` or `_` boundary immediately after `filename`
    so that a prefix like `1_S1` does not match `1_S10*` / `1_S100*` etc.
    Without the anchor, two samples sharing a numeric prefix would silently
    cross-pair (or trip the "Multiple R1 files found" error path).
    """
    data_dir = Path(data_dir).resolve()
    R1, R2 = None, None

    for file in data_dir.glob(f"{filename}[._]*"):
        name = file.name
        if _FASTQ_R1_PATTERN.search(name):
            if R1 is not None:
                raise ValueError(
                    f"Multiple R1 files found for prefix '{filename}': {R1.name} and {name}"
                )
            R1 = file
        elif _FASTQ_R2_PATTERN.search(name):
            if R2 is not None:
                raise ValueError(
                    f"Multiple R2 files found for prefix '{filename}': {R2.name} and {name}"
                )
            R2 = file

    if not R1 or not R2:
        raise FileNotFoundError(
            f"Could not find matching R1 and R2 fastq files for prefix '{filename}' in {data_dir}. "
            f"Expected files matching patterns like {filename}_R1_001.fastq.gz, {filename}_R1.fq.gz, "
            f"{filename}_1.fastq, etc."
        )

    return R1, R2


def file_digest(path, length=12):
    """Return a short hex digest of a file's contents.

    Used to key the bbmap index directory by reference content, so that a
    changed reference produces a different output path and bbmap is never
    pointed at an index built from a different reference. The risk this
    avoids is silent: bbmap with `rebuild=t` rewrites its own index files
    but doesn't necessarily clean up unrelated artifacts from a previous
    build, which can confuse mapping.

    SHA-1 is fine here — we're keying a filesystem path, not authenticating.
    """
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()[:length]


def run_script(snakemake, func):
    """Set up logging from Snakemake log config, run func, and log any unhandled exceptions."""
    log_file = snakemake.log[0]
    if log_file:
        logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG)
    try:
        func(snakemake)
    except Exception:
        logging.exception("Script failed")
        raise


def set_index_with_unique_check(
    dataframe: pd.DataFrame, index_column: str, drop: bool = True
) -> pd.DataFrame:
    """Set a DataFrame index after enforcing uniqueness in the source column.
    Lists duplicate values in case of error."""
    duplicated_index = dataframe[index_column].duplicated(keep=False)

    if duplicated_index.any():
        duplicate_values = sorted(
            dataframe.loc[duplicated_index, index_column].unique()
        )
        raise ValueError(
            f"Duplicate values found in '{index_column}': {duplicate_values}"
        )
    return dataframe.set_index(index_column, drop=drop)
