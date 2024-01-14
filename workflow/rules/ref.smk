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
    """Index the reference sequence using samtools. Necessary for GATK.
    
    This is a non-wrapper version of the rule to allow local installs of samtools."""
    input:
        expand("{ref_dir}/{{reference}}.fasta", ref_dir=config["ref_dir"]),
    output:
        expand("{ref_dir}/{{reference}}.fasta.fai", ref_dir=config["ref_dir"]),
    log:
        "logs/gatk/{reference}.create_index.log",
    shell:
        "samtools faidx {input} 2> {log}"

rule prepare_bbmap_index:
    """Generate the index for mapping with bbmap. This must be run once before mapping."""
    input:
        expand(
            "{ref_dir}/{reference}",
            ref_dir=config["ref_dir"],
            reference=config["reference"],
        ),
    output:
        "ref/genome/1/chr1.chrom.gz",
    threads: 16
    log:
        expand(
            "logs/{experiment}/bbmap/{reference}.bbmap_index.log",
            experiment=experiment,
            reference=config["reference"],
        ),
    shell:
        "bbmap.sh "
        "ref={input}"
