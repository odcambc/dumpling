rule trim_clean_correct:
    """Stream the three BBTools read-processing stages (adapter trim, contaminant
    removal, error correction) through one pipeline so the two transient FASTQ
    pairs between them never hit disk. Only the final `.ec.clean.trim.fastq.gz`
    pair is written. Each stage still writes its own stats files in parallel
    with the streamed payload, so MultiQC's BBDuk / BBMerge sections and the
    baseline file list remain unchanged.

    bbmap_map (in map.smk) stays a separate rule because its 8.5GB RSS would
    push the merged rule into double-digit per-sample memory and `_map.covstats`
    is the only stats file that's a load-bearing rule input downstream."""
    input:
        unpack(get_file_from_sample),
    output:
        R1_ec=temp("results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz"),
        R2_ec=temp("results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz"),
        # Per-tool stats files preserved at their historical paths so MultiQC's
        # auto-discovery and generate_baseline_file_list.py both keep working.
        qhist="stats/{experiment}/{sample_prefix}_trim.qhist",
        bhist="stats/{experiment}/{sample_prefix}_trim.bhist",
        gchist="stats/{experiment}/{sample_prefix}_trim.gchist",
        aqhist="stats/{experiment}/{sample_prefix}_trim.aqhist",
        lhist="stats/{experiment}/{sample_prefix}_trim.lhist",
        trim_stats="stats/{experiment}/{sample_prefix}_trim.stats.txt",
        contam_stats="stats/{experiment}/{sample_prefix}_trim_contam.stats.txt",
        ihist="stats/{experiment}/{sample_prefix}_merge.ihist",
    params:
        adapters=adapters_ref,
        contaminants=contaminants_ref,
        mem=config["mem"],
        compression_flags=bbtools_compression_flags,
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.trim_clean_correct.benchmark.txt"
    log:
        bbduk_trim="logs/{experiment}/bbduk/{sample_prefix}.trim.bbduk.log",
        bbduk_clean="logs/{experiment}/bbduk/{sample_prefix}.clean.bbduk.log",
        bbmerge="logs/{experiment}/bbmerge/{sample_prefix}.bbmerge.log",
    threads: 16
    shell:
        # Three BBTools JVMs run concurrently in a pipe. Each is given a 5/5/6
        # thread budget so the total stays at `threads`; passing t={threads} to
        # all three would oversubscribe (3*16=48 software threads on 16 cores).
        # Inter-stage payload is uncompressed interleaved FASTQ on stdin/stdout;
        # only the final bbmerge write to disk uses the compression flags.
        "( bbduk.sh "
        "-Xmx{params.mem}g "
        "in1={input.R1} in2={input.R2} "
        "ref={params.adapters} ktrim=r k=23 mink=10 hdist=1 tpe tbo "
        "out=stdout.fq interleaved=t "
        "bhist={output.bhist} "
        "qhist={output.qhist} "
        "gchist={output.gchist} "
        "aqhist={output.aqhist} "
        "lhist={output.lhist} "
        "stats={output.trim_stats} "
        "overwrite=false t=5 gcbins=auto 2> {log.bbduk_trim} "
        "| bbduk.sh "
        "-Xmx{params.mem}g "
        "in=stdin.fq interleaved=t "
        "ref={params.contaminants} k=31 "
        "out=stdout.fq interleaved=t "
        "stats={output.contam_stats} "
        "overwrite=false t=5 2> {log.bbduk_clean} "
        "| bbmerge.sh "
        "-Xmx{params.mem}g "
        "in=stdin.fq interleaved=t "
        "out1={output.R1_ec} out2={output.R2_ec} "
        "ihist={output.ihist} "
        "ecco mix showhiststats=t "
        "{params.compression_flags}"
        "overwrite=false t=6 2> {log.bbmerge} "
        ")"
