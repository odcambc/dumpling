"""Load + validate the barcode->variant map once and write the duplicates report.

Ambiguous tags (one barcode naming multiple distinct variants) are dropped by
load_barcode_map's drop-all policy and recorded here. This runs as a single job
with no sample wildcard, so the parallel per-sample count_barcodes jobs never
race on the report file. The declared output is always written — empty (header
only) when the map is clean — so Snakemake sees the output exist.
"""

import logging
from pathlib import Path

from script_utils import load_barcode_map, run_script


def _run(snakemake):
    barcode_map_file = snakemake.config["barcode_map"]
    report = Path(snakemake.output.duplicates)

    barcode_to_variant, conflicts = load_barcode_map(barcode_map_file)

    report.parent.mkdir(parents=True, exist_ok=True)
    conflicts.to_csv(report, index=False)

    logging.info(
        "Barcode map %s: %d tags kept, %d ambiguous tag(s) dropped -> %s",
        barcode_map_file,
        len(barcode_to_variant),
        len(conflicts),
        report,
    )


def main():
    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
