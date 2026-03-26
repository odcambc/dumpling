import logging

import pandas as pd


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
