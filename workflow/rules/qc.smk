def return_directories(wildcards):
    """This function takes a list of files and returns the unique directories the are in."""
    directories = [
        "stats/" + wildcards.experiment_name + "/fastqc/",
        "results/" + wildcards.experiment_name + "/gatk/",
        "stats/" + wildcards.experiment_name + "/",
    ]

    return directories


rule multiqc_dir:
    """Final QC: aggregate FastQC and intermediate log files into a final report with MultiQC.
  """
    input:
        return_directories,
    output:
        "stats/{experiment_name}/{experiment_name}_multiqc.html",
    benchmark:
        "benchmarks/{experiment_name}/multiqc.benchmark.txt"
    shell:
        "multiqc "
        "--outdir stats/{wildcards.experiment_name}/ "
        "--filename {wildcards.experiment_name}_multiqc.html "
        "{input}"


rule fastqc:
    """Initial QC: run FastQC on all input reads.
  """
    input:
        expand("{data_dir}/{{sample_read}}_001.fastq.gz", data_dir=config["data_dir"]),
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
