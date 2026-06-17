# The map/GATK-stage stats that multiqc_dir depends on only exist in direct
# (align -> GATK) mode. Barcode mode bypasses both stages, so those triggers are
# dropped and MultiQC aggregates only the read-level + barcode-counting stats
# (FastQC, trim/clean/correct, total_processing). In direct mode the map-stage
# trigger depends on the aligner: BBMap emits `_map.covstats`; minimap2 emits
# `_samtools_stats.txt`. MultiQC autodiscovers everything else by content.
if config["barcoded"]:
    _stage_stats_triggers = []
else:
    if config["aligner"] == "bbmap":
        _map_stats_trigger = "stats/{{experiment}}/{sample_prefix}_map.covstats"
    else:  # minimap2
        _map_stats_trigger = "stats/{{experiment}}/{sample_prefix}_samtools_stats.txt"
    _stage_stats_triggers = expand(
        _map_stats_trigger, sample_prefix=samples
    ) + expand(
        "results/{{experiment}}/gatk/{sample_prefix}.variantCounts",
        sample_prefix=samples,
    )


rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC."""
    input:
        _stage_stats_triggers,
        [
            f"stats/{{experiment}}/fastqc/{fastqc_names[f]['R1']}_fastqc.html"
            for f in files
        ],
        [
            f"stats/{{experiment}}/fastqc/{fastqc_names[f]['R2']}_fastqc.html"
            for f in files
        ],
        expand(
            "stats/{{experiment}}/processing/{sample_prefix}_total_processing.tsv",
            sample_prefix=samples,
        ),
    output:
        "stats/{experiment}/{experiment}_multiqc.html",
    log:
        "logs/{experiment}/multiqc.log",
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    params:
        extra="-c config/multiqc_config.yaml",
    wrapper:
        "v3.1.0/bio/multiqc"


rule fastqc:
    """Initial QC: run FastQC on all input reads."""
    input:
        lambda wc: fastqc_input_map[wc.fastqc_name],
    output:
        html="stats/{experiment}/fastqc/{fastqc_name}_fastqc.html",
        zip="stats/{experiment}/fastqc/{fastqc_name}_fastqc.zip",
    log:
        "logs/{experiment}/fastqc/{fastqc_name}.log",
    benchmark:
        "benchmarks/{experiment}/{fastqc_name}.fastqc.benchmark.txt"
    threads: 8
    resources:
        mem_mb=config["mem_fastqc"],
    params:
        "--quiet",
    wrapper:
        "v3.1.0/bio/fastqc"
