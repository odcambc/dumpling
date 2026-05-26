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


if config["aligner"] == "bbmap":

    rule prepare_bbmap_index:
        """Generate the index for mapping with bbmap. This must be run once before mapping.

        The output directory is keyed on a hash of the reference contents
        (`{experiment}_{ref_digest}`), so a changed reference produces a fresh
        directory and bbmap is never pointed at an index built from a different
        reference. No pre-build cleanup is needed: a different reference goes
        somewhere else by construction, and a re-run against the same reference
        rebuilds via `rebuild=t`."""
        input:
            str(reference_file),
        output:
            index_dir=temp(directory(f"ref/bbmap/{experiment}_{ref_digest}")),
        params:
            mem=config["mem"],
            kmers=config["kmers"],
            compression_flags=bbtools_compression_flags,
        threads: 16
        log:
            f"logs/{experiment}/bbmap/{reference_name}.bbmap_index.log",
        shell:
            "bbmap.sh "
            "-Xmx{params.mem}g "
            "ref={input} "
            "path={output.index_dir} "
            "build=1 "
            "{params.compression_flags}"
            "k={params.kmers} "
            "rebuild=t "
            "2> {log}"

elif config["aligner"] == "minimap2":

    rule prepare_minimap2_index:
        """Pre-build the minimap2 short-read index (`-x sr`) so each per-sample
        mapping doesn't re-index the reference. Index is keyed on the ref
        digest so a changed reference produces a fresh file by construction."""
        input:
            str(reference_file),
        output:
            index=temp(f"ref/minimap2/{experiment}_{ref_digest}.mmi"),
        threads: 4
        log:
            f"logs/{experiment}/minimap2/{reference_name}.minimap2_index.log",
        shell:
            "mkdir -p $(dirname {output.index}); "
            "minimap2 -x sr -t {threads} -d {output.index} {input} 2> {log}"
