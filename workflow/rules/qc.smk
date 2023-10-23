rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC.
  """
    input:
        expand("stats/{{experiment_name}}/{sample_prefix}_map.covstats",
        sample_prefix=samples),
        expand("stats/{{experiment_name}}/fastqc/{file}_R1_001_fastqc.html",
        file=files),
        expand("results/{{experiment_name}}/gatk/{sample_prefix}.variantCounts",
        sample_prefix=samples)
    output:
        "stats/{experiment_name}/{experiment_name}_multiqc.html",
    benchmark:
        "benchmarks/{experiment_name}/multiqc.benchmark.txt"
    log:
        "logs/{experiment_name}/multiqc.log"
    wrapper:
        "v2.6.0/bio/multiqc"


rule fastqc:
    """Initial QC: run FastQC on all input reads.
  """
    input:
        expand("{data_dir}/{{sample_read}}.fastq.gz", data_dir=config["data_dir"]),
    output:
        html="stats/{experiment}/fastqc/{sample_read}_fastqc.html",
        zip="stats/{experiment}/fastqc/{sample_read}_fastqc.zip",
    params:
        "--quiet",
    benchmark:
        "benchmarks/{experiment}/{sample_read}.fastqc.benchmark.txt"
    log:
        "logs/{experiment}/fastqc/{sample_read}.log",
    wrapper:
        "v2.6.0/bio/fastqc"
