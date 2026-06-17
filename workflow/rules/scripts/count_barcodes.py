"""Count barcode reads into the per-sample enrich-format TSV (barcode mode).

Barcode-counting mode replaces the align -> GATK -> process_sample path: rather
than calling a mutation from aligned reads, each read carries a short tag that a
user-provided map associates with a variant. We extract the tag from R1 (already
overlap-error-corrected against R2 by trim_clean_correct's `bbmerge ... ecco`
step), match it exactly against the map whitelist, and tally it. The outputs are
the same `(hgvs, count)` enrich-format TSV + processed CSV + total-stats file
that process_sample emits in direct mode, so every downstream scoring/deposit
rule is reused unchanged.

WT/synonymous handling: the enrich writer folds synonymous-variant counts into
the `_wt` pseudovariant the scoring backends normalize against. Barcode mode has
no GATK-derived `mutation_type`, so we recognize synonymous variants from their
HGVS (the same convention format_mavedb uses: `p.(=)`, `p.(L54=)`, or same-AA
`p.(L54L)`) and synthesize the `mutation_type` column write_enrich_df keys off.
"""

import gzip
import logging
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import process_variants
from script_utils import load_barcode_map, run_script

# Mirrors format_mavedb._SYNONYMOUS_RE: a same-AA substitution p.(L54L) or the
# explicit p.(L54=) / population-level p.(=) form. Used to fold synonymous
# counts into the _wt pseudovariant exactly as write_enrich_df does for
# mutation_type == "S".
_SYNONYMOUS_RE = re.compile(r"^p\.\(([A-Z*])(\d+)([A-Z*=])\)$")


def is_synonymous(hgvs):
    """True if an HGVS protein string denotes no amino-acid change."""
    if hgvs == "p.(=)":
        return True
    m = _SYNONYMOUS_RE.match(hgvs)
    return bool(m) and (m.group(3) == "=" or m.group(1) == m.group(3))


def make_extractor(barcode_pattern=None, barcode_start=None, barcode_length=None):
    """Return a ``read_seq -> barcode|None`` extractor for the configured spec.

    Exactly one of (regex pattern) / (start+length) is set — validate_barcode_config
    enforces that upstream. Returns None when the tag can't be pulled: a regex
    miss, or a read shorter than the fixed-position window.
    """
    if barcode_pattern:
        compiled = re.compile(barcode_pattern)

        def extract(seq):
            match = compiled.search(seq)
            return match.group(1) if match else None

        return extract

    start = barcode_start
    end = barcode_start + barcode_length

    def extract(seq):
        if len(seq) < end:
            return None
        return seq[start:end]

    return extract


def iter_fastq_sequences(fastq_path):
    """Yield uppercased sequence lines from a (optionally gzipped) FASTQ."""
    path = str(fastq_path)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:  # FASTQ: the sequence is the 2nd of each 4-line record
                yield line.strip().upper()


def count_barcodes(fastq_path, barcode_to_variant, extract, min_barcode_count=0):
    """Stream reads, extract + exact-match tags, return (barcode_counts, stats).

    stats is a Counter of read-disposition tallies (total / no-barcode / unknown
    / matched / below-min) for the per-sample QC report.
    """
    barcode_counts = Counter()
    stats = Counter()

    for seq in iter_fastq_sequences(fastq_path):
        stats["total_reads"] += 1
        barcode = extract(seq)
        if barcode is None:
            stats["no_barcode_extracted"] += 1
            continue
        if barcode not in barcode_to_variant:
            stats["unknown_barcode"] += 1
            continue
        barcode_counts[barcode] += 1

    stats["matched_reads"] = sum(barcode_counts.values())

    if min_barcode_count > 0:
        below = [bc for bc, c in barcode_counts.items() if c < min_barcode_count]
        stats["below_min_count_barcodes"] = len(below)
        stats["below_min_count_reads"] = sum(barcode_counts[bc] for bc in below)
        for bc in below:
            del barcode_counts[bc]

    stats["distinct_barcodes_observed"] = len(barcode_counts)
    return barcode_counts, stats


def aggregate_variant_counts(barcode_counts, barcode_to_variant):
    """Sum surviving per-barcode counts into per-variant counts over the FULL
    map variant universe (0 for unobserved), so every sample shares one variant
    set — scoring backends compare counts across samples and need a stable set.
    """
    variant_counts = {variant: 0 for variant in barcode_to_variant.values()}
    for barcode, count in barcode_counts.items():
        variant_counts[barcode_to_variant[barcode]] += count
    return variant_counts


def build_counts_frame(variant_counts):
    """Build the (hgvs, count, mutation_type) frame write_enrich_df consumes.

    The synthesized mutation_type ("S" for synonymous HGVS, else "M") is what
    write_enrich_df uses to compute the _wt synonymous-sum row, so the emitted
    enrich TSV matches the direct-mode shape.
    """
    df = pd.DataFrame(
        {"hgvs": list(variant_counts), "count": list(variant_counts.values())}
    )
    df["mutation_type"] = df["hgvs"].apply(lambda h: "S" if is_synonymous(h) else "M")
    return df


def _run(snakemake):
    config = snakemake.config
    sample_name = snakemake.wildcards.sample_prefix
    fastq_path = snakemake.input.r1
    barcode_map_file = config["barcode_map"]
    min_barcode_count = config.get("min_barcode_count", 0)

    # Reload the (already validated) map per job — the drop-all conflict policy
    # is deterministic, so every sample gets the same clean whitelist; the
    # duplicates report is written once by the barcode_map_report rule.
    barcode_to_variant, _conflicts = load_barcode_map(barcode_map_file)
    logging.info(
        "Loaded %d barcodes -> %d distinct variants from %s",
        len(barcode_to_variant),
        len(set(barcode_to_variant.values())),
        barcode_map_file,
    )

    extract = make_extractor(
        barcode_pattern=config.get("barcode_pattern"),
        barcode_start=config.get("barcode_start"),
        barcode_length=config.get("barcode_length"),
    )

    barcode_counts, stats = count_barcodes(
        fastq_path, barcode_to_variant, extract, min_barcode_count
    )
    logging.info(
        "Sample %s: %d/%d reads matched the whitelist (%d unknown, %d no-tag)",
        sample_name,
        stats["matched_reads"],
        stats["total_reads"],
        stats["unknown_barcode"],
        stats["no_barcode_extracted"],
    )

    variant_counts = aggregate_variant_counts(barcode_counts, barcode_to_variant)
    counts_df = build_counts_frame(variant_counts)

    # enrich-format TSV (the seam) + processed CSV (mavedb counts: hgvs, count).
    process_variants.write_enrich_df(snakemake.output.enrich, counts_df, noprocess=False)

    csv_path = Path(snakemake.output.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    counts_df[["hgvs", "count"]].to_csv(csv_path, index=False)

    stats["variants_in_map"] = len(variant_counts)
    process_variants.write_stats_file(snakemake.output.total_stats, dict(stats))


def main():
    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
