import hashlib
import logging
from pathlib import Path

import pandas as pd


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
