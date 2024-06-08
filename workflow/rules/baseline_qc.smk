rule multiqc_baseline:
    """Rule for generating baseline QC report with MultiQC.
    """
    input:
        expand("results/{{experiment}}/processed_counts/{samples}.tsv", samples=samples),
        expand(
            "stats/{{experiment}}/processing/{samples}_total_processing.tsv",
            samples=samples,
        ),
        file_list="stats/{experiment}/baseline_file_list.txt",
        config="config/{experiment}_multiqc_baseline_config.yaml",
    output:
        "stats/{experiment}/{experiment}_baseline_multiqc.html",
    benchmark:
        "benchmarks/{experiment}/multiqc.baseline.benchmark.txt"
    params:
        outdir="stats/{experiment}/",
        outfile="{experiment}_baseline_multiqc.html",
    log:
        "logs/{experiment}/multiqc.baseline.log",
    conda:
        "../envs/multiqc-baseline.yaml"
    shell:
        "multiqc "
        "-c {input.config} "
        "--file-list {input.file_list} "
        "--outdir {params.outdir} "
        "--filename {params.outfile} "
        "--force "
        "--verbose "
        "--interactive "
        "2>{log}"


rule generate_baseline_configs:
    """Rule for generating a multiQC config file for baseline QC."""
    output:
        "config/{experiment}_multiqc_baseline_config.yaml",
    log:
        "logs/{experiment}/multiqc_baseline_config.log",
    script:
        "scripts/generate_baseline_configs.py"


rule generate_baseline_file_list:
    """Rule for generating a file list for baseline QC."""
    input:
        expand(
            "stats/{{experiment}}/fastqc/{file}_R1_001_fastqc.html",
            file=baseline_files,
        ),
        expand(
            "stats/{{experiment}}/fastqc/{file}_R2_001_fastqc.html",
            file=baseline_files,
        ),
        expand(
            "results/{{experiment}}/gatk/{sample_prefix}.variantCounts",
            sample_prefix=baseline_samples,
        ),
    output:
        "stats/{experiment}/baseline_file_list.txt",
    log:
        "logs/{experiment}/multiqc_baseline_file_list.log",
    script:
        "scripts/generate_baseline_file_list.py"
