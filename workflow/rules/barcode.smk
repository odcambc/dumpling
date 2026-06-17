# Barcode-counting mode (config["barcoded"]). These rules replace the
# align -> GATK -> process_sample path with a direct tag-extraction-and-tally
# step, converging on the same per-sample enrich-format TSV. They are defined
# only in barcode mode; process_sample (process.smk) is gated off in the same
# mode so exactly one rule produces the shared enrich_format/{sample}.tsv path.
if config["barcoded"]:

    rule barcode_map_report:
        """Load + validate the barcode->variant map once and emit the duplicates
        report for any ambiguous (multi-variant) tags dropped by the drop-all
        policy. One job (no sample wildcard) so the per-sample count_barcodes
        jobs don't race on the report file."""
        output:
            duplicates="results/{experiment}/duped_barcodes.csv",
        log:
            "logs/{experiment}/barcode_map_report.log",
        script:
            "scripts/barcode_map_report.py"

    rule count_barcodes:
        """Barcode-mode producer of the per-sample enrich-format counts: extract
        the tag from each (already overlap-error-corrected) R1 read, exact-match
        it against the map whitelist, and tally into the same (hgvs, count) TSV +
        processed CSV + total-stats file process_sample emits in direct mode.
        Mutually exclusive with process_sample (both gate on config['barcoded'])."""
        input:
            r1="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
            # Declared so edits to the map re-trigger counting. Only referenced
            # in barcode mode, where validate_barcode_config guarantees it exists.
            barcode_map=config["barcode_map"],
        output:
            enrich="results/{experiment}/processed_counts/enrich_format/{sample_prefix}.tsv",
            csv="results/{experiment}/processed_counts/{sample_prefix}.csv",
            total_stats="stats/{experiment}/processing/{sample_prefix}_total_processing.tsv",
        resources:
            mem_mb=config["mem_process_sample"],
        benchmark:
            "benchmarks/{experiment}/{sample_prefix}.count_barcodes.benchmark.txt"
        log:
            "logs/{experiment}/scripts/{sample_prefix}.count_barcodes.log",
        script:
            "scripts/count_barcodes.py"
