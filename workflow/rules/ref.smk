rule prepare_reference:
    """Normalize the reference FASTA: copy the user's reference into a working
    path under `ref/` and rewrite the first `>` line to `>{reference_name}` so
    the contig name inside the FASTA matches the filename stem. Solves
    Upstream #19, where reference generation fails when the user's FASTA
    header doesn't match `Path(config["reference"]).stem` (BBMap index build,
    GATK SAM/ref contig matching, etc.).

    Downstream index/dict/aligner rules all consume from this normalized
    path; the user's original file is left untouched."""
    input:
        str(reference_file),
    output:
        str(normalized_reference_file),
    log:
        f"logs/{reference_name}.prepare_reference.log",
    shell:
        # awk rewrites only the FIRST `>` line; subsequent contigs are passed
        # through unchanged. DMS references are single-contig in practice, but
        # leaving multi-contig FASTAs intact avoids accidentally renaming
        # secondary contigs the user might be relying on. Tee'ing the result
        # to the log file gives a copy of the rewritten contig name for
        # debugging.
        "mkdir -p $(dirname {output}); "
        "awk -v name='{reference_name}' "
        "    'NR==1 && /^>/ {{ sub(/^>.*/, \">\" name); }} {{ print }}' "
        "    {input} > {output} 2> {log}; "
        "echo \"Normalized $(head -1 {input}) -> $(head -1 {output})\" >> {log}"


rule prepare_dict_GATK:
    """Generate a sequence dictionary from the reference sequence. Necessary for GATK.

    https://gatk.broadinstitute.org/hc/en-us/articles/360036712531-CreateSequenceDictionary-Picard-
    """
    input:
        str(normalized_reference_file),
    output:
        str(normalized_reference_file.parent / f"{reference_name}.dict"),
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
            str(normalized_reference_file),
        output:
            str(normalized_reference_file) + ".fai",
        log:
            f"logs/gatk/{reference_name}.create_index.log",
        shell:
            "samtools faidx {input} 2> {log}"

else:

    rule prepare_index:
        """Index the reference sequence using samtools. Necessary for GATK."""
        input:
            str(normalized_reference_file),
        output:
            str(normalized_reference_file) + ".fai",
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
            str(normalized_reference_file),
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
            str(normalized_reference_file),
        output:
            index=temp(f"ref/minimap2/{experiment}_{ref_digest}.mmi"),
        threads: 4
        log:
            f"logs/{experiment}/minimap2/{reference_name}.minimap2_index.log",
        shell:
            "mkdir -p $(dirname {output.index}); "
            "minimap2 -x sr -t {threads} -d {output.index} {input} 2> {log}"
