rule multiqc_dir:
    input:
        "stats/fastqc/{sample}_R1_fastqc.html"
    output:
        "stats/multiqc.html"
    params:
        extra=""
    benchmark:
        "benchmarks/multiqc.benchmark.txt"
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"

rule fastqc:
  input:
    expand("{data_dir}/{experiment_dir}/{{sample}}_001.fastq.gz", data_dir=config["data_dir"], experiment_dir=config["experiment_dir"])
  output:
    html="stats/fastqc/{sample}_fastqc.html",
    zip="stats/fastqc/{sample}_fastqc.zip"
  params: "--quiet"
  benchmark:
    "benchmarks/{sample}.fastqc.benchmark.txt"
  log:
    "logs/fastqc/{sample}.log"
  wrapper:
    "v1.5.0/bio/fastqc"
