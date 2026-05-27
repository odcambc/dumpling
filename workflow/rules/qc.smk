# The map-stage stats file that triggers the multiqc_dir dependency on the
# mapping rule depends on which aligner we're using. BBMap emits
# `_map.covstats`; minimap2's rule emits `_samtools_stats.txt` (via
# `samtools stats` on the BAM). MultiQC autodiscovers everything else in
# the stats directory by content.
if config["aligner"] == "bbmap":
    _map_stats_trigger = "stats/{{experiment}}/{sample_prefix}_map.covstats"
else:  # minimap2
    _map_stats_trigger = "stats/{{experiment}}/{sample_prefix}_samtools_stats.txt"


rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC.
  """
    input:
        expand(_map_stats_trigger, sample_prefix=samples),
        [
            f"stats/{{experiment}}/fastqc/{fastqc_names[f]['R1']}_fastqc.html"
            for f in files
        ],
        [
            f"stats/{{experiment}}/fastqc/{fastqc_names[f]['R2']}_fastqc.html"
            for f in files
        ],
        expand(
            "results/{{experiment}}/gatk/{sample_prefix}.variantCounts",
            sample_prefix=samples,
        ),
        expand(
            "stats/{{experiment}}/processing/{sample_prefix}_total_processing.tsv",
            sample_prefix=samples,
        ),
    output:
        "stats/{experiment}/{experiment}_multiqc.html",
    params:
        extra="-c config/multiqc_config.yaml",
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    log:
        "logs/{experiment}/multiqc.log",
    wrapper:
        "v3.1.0/bio/multiqc"


rule fastqc:
    """Initial QC: run FastQC on all input reads."""
    input:
        lambda wc: fastqc_input_map[wc.fastqc_name],
    output:
        html="stats/{experiment}/fastqc/{fastqc_name}_fastqc.html",
        zip="stats/{experiment}/fastqc/{fastqc_name}_fastqc.zip",
    params:
        "--quiet",
    benchmark:
        "benchmarks/{experiment}/{fastqc_name}.fastqc.benchmark.txt"
    log:
        "logs/{experiment}/fastqc/{fastqc_name}.log",
    threads: 8
    resources:
        mem_mb=config["mem_fastqc"],
    wrapper:
        "v3.1.0/bio/fastqc"
