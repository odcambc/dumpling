rule multiqc_dir:
  """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC.
  """

  input:
      "stats/{experiment}/fastqc/",
      "stats/{experiment}/",
      "results/{experiment}/gatk/"
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
  """Initial QC: run FastQC on all input reads.
  """
  input:
    expand("{data_dir}/{{sample_read}}_001.fastq.gz", data_dir=config["data_dir"])
  output:
    html="stats/{experiment}/fastqc/{sample_read}_fastqc.html",
    zip="stats/{experiment}/fastqc/{sample_read}_fastqc.zip"
  params: "--quiet"
  benchmark:
    "benchmarks/{experiment}/{sample_read}.fastqc.benchmark.txt"
  log:
    "logs/{experiment}/fastqc/{sample_read}.log"
  wrapper:
    "v1.5.0/bio/fastqc"