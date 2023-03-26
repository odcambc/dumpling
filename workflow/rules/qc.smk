rule multiqc_dir:
    input:
        "stats/{experiment}/fastqc/{sample}_R1_fastqc.html"
    output:
        "stats/{experiment}/multiqc.html"
    params:
        extra=""
    benchmark:
        "benchmarks/{experiment}/multiqc.benchmark.txt"
    log:
        "logs/{experiment}/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"

rule fastqc:
  input:
    expand("{data_dir}/{{sample_name}}_001.fastq.gz", data_dir=config["data_dir"])
  output:
    html="stats/{experiment}/fastqc/{sample_name}_fastqc.html",
    zip="stats/{experiment}/fastqc/{sample_name}_fastqc.zip"
  params: "--quiet"
  benchmark:
    "benchmarks/{experiment}/{sample_name}.fastqc.benchmark.txt"
  log:
    "logs/{experiment}/fastqc/{sample_name}.log"
  wrapper:
    "v1.5.0/bio/fastqc"