rule prepare_dict_GATK:
    """Generate a sequence dictionary from the reference sequence. Necessary for GATK.

    https://gatk.broadinstitute.org/hc/en-us/articles/360036712531-CreateSequenceDictionary-Picard-
    """
    input:
        expand("{ref_dir}/{{reference}}.fasta", ref_dir=config["ref_dir"]),
    output:
        expand("{ref_dir}/{{reference}}.dict", ref_dir=config["ref_dir"]),
    log:
        "logs/gatk/{reference}.create_dict.log",
    shell:
        "gatk CreateSequenceDictionary -R {input} 2> {log}"


rule prepare_index:
    """Index the reference sequence using samtools. Necessary for GATK."""
    input:
        expand("{ref_dir}/{{reference}}.fasta", ref_dir=config["ref_dir"]),
    output:
        expand("{ref_dir}/{{reference}}.fasta.fai", ref_dir=config["ref_dir"]),
    log:
        "logs/gatk/{reference}.create_index.log",
    wrapper:
        "v2.8.0/bio/samtools/faidx"