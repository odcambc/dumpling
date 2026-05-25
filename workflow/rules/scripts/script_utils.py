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


def load_experiments(experiment_file) -> pd.DataFrame:
    """Load the experiment-definition CSV with the canonical encoding and
    indexing the rest of the pipeline expects.

    Every script that needs to look up samples calls this so the read
    semantics (utf-8-sig BOM tolerance, blank-row drop, unique-sample
    enforcement) stay identical across rules. Drift between callers
    previously meant a CSV that worked for one script could break another.
    """
    df = pd.read_csv(experiment_file, header=0, encoding="utf-8-sig").dropna(
        how="all"
    )
    return set_index_with_unique_check(df, "sample", drop=False)


def translate_orf(ref_sequence, orf_range: str):
    """Translate the ORF region of a reference DNA sequence.

    The ``orf_range`` string follows ``"<start>-<end>"`` with 1-based,
    inclusive positions. Two conventions are supported:

    - **Linear** (``start <= end``): returns ``ref_sequence[start-1:end]`` translated.
    - **Circular** (``start > end``, i.e. the ORF crosses the origin):
      returns ``ref_sequence[start-1:] + ref_sequence[:end]`` translated.

    Lives in script_utils so generate_variants.py and process_counts.py
    produce the same amino-acid reference from the same config — they
    used to diverge for circular ORFs, with process_counts silently
    returning an empty translation while generate_variants did the
    correct concatenation.
    """
    orf_start, orf_end = map(int, orf_range.split("-"))
    if orf_start > orf_end:
        nt = ref_sequence[orf_start - 1 :] + ref_sequence[:orf_end]
    else:
        nt = ref_sequence[orf_start - 1 : orf_end]
    return nt.translate()


def validate_experiment_time_or_bin(experiments: pd.DataFrame) -> None:
    """Each experiment row must set exactly one of ``time`` (timecourse) or
    ``bin`` (FACS-type). The JSON-Schema ``oneOf`` enforcing this raised
    "is not valid under any of the given schemas" without naming the offending
    sample, which was effectively unactionable for users with a 100-row CSV.
    """
    has_time_col = "time" in experiments.columns
    has_bin_col = "bin" in experiments.columns

    problems = []
    for idx, row in experiments.iterrows():
        time_set = has_time_col and pd.notnull(row.get("time"))
        bin_set = has_bin_col and pd.notnull(row.get("bin"))
        sample = row.get("sample", f"<row {idx}>")
        if not time_set and not bin_set:
            problems.append(
                f"sample {sample!r}: missing 'time' or 'bin' — set 'time' "
                f"for timecourse experiments or 'bin' for FACS-type experiments"
            )
        elif time_set and bin_set:
            problems.append(
                f"sample {sample!r}: 'time' and 'bin' are both set — pick "
                f"one (timecourse uses 'time', FACS-type uses 'bin')"
            )

    if problems:
        details = "\n  - ".join(problems)
        raise ValueError(
            "Invalid experiment definitions:\n  - " + details
        )
