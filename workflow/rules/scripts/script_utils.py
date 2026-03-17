import logging


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
