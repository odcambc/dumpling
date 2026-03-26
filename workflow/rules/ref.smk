rule prepare_dict_GATK:
    """Generate a sequence dictionary from the reference sequence. Necessary for GATK.

    https://gatk.broadinstitute.org/hc/en-us/articles/360036712531-CreateSequenceDictionary-Picard-
    """
    input:
        str(reference_file),
    output:
        str(reference_file.parent / f"{reference_name}.dict"),
    log:
        f"logs/gatk/{reference_name}.create_dict.log",
    shell:
        "gatk CreateSequenceDictionary -R {input} 2> {log}"


if samtools_local:

    rule prepare_index_nowrapper:
        """Index the reference sequence using samtools. Necessary for GATK.

        This uses the samtools in the path rather than the wrapper.
        This can be useful if conda does not have a samtools version for the platform.
        """
        input:
            str(reference_file),
        output:
            str(reference_file) + ".fai",
        log:
            f"logs/gatk/{reference_name}.create_index.log",
        shell:
            "samtools faidx {input} 2> {log}"

else:

    rule prepare_index:
        """Index the reference sequence using samtools. Necessary for GATK."""
        input:
            str(reference_file),
        output:
            str(reference_file) + ".fai",
        log:
            f"logs/gatk/{reference_name}.create_index.log",
        wrapper:
            "v3.12.0/bio/samtools/faidx"


rule prepare_bbmap_index:
    """Generate the index for mapping with bbmap. This must be run once before mapping."""
    input:
        str(reference_file),
    output:
        index_dir=temp(directory(f"ref/bbmap/{experiment}")),
    params:
        mem=config["mem"],
        kmers=config["kmers"],
        compression_flags=bbtools_compression_flags,
    threads: 16
    log:
        f"logs/{experiment}/bbmap/{reference_name}.bbmap_index.log",
    shell:
        "rm -rf {output.index_dir} "
        "&& mkdir -p {output.index_dir} "
        "&& "
        "bbmap.sh "
        "-Xmx{params.mem}g "
        "ref={input} "
        "path={output.index_dir} "
        "build=1 "
        "{params.compression_flags}"
        "k={params.kmers} "
        "rebuild=t "
        "2> {log}"
