rule map_to_reference_bbmap:
    """Map reads to reference sequence using BBMap, streaming the SAM through
    `samtools view -b` so the on-disk intermediate is a BAM. covhist output is
    currently disabled as it causes MultiQC bloat."""
    input:
        index_dir=f"ref/bbmap/{experiment}_{ref_digest}",
        R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
        R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
    output:
        bam=temp("results/{experiment}/{sample_prefix}.mapped.bam"),
        covstats="stats/{experiment}/{sample_prefix}_map.covstats",
        basecov="stats/{experiment}/{sample_prefix}_map.basecov",
        bincov="stats/{experiment}/{sample_prefix}_map.bincov",
        ehist="stats/{experiment}/{sample_prefix}_map.ehist",
        indelhist="stats/{experiment}/{sample_prefix}_map.indelhist",
        mhist="stats/{experiment}/{sample_prefix}_map.mhist",
        idhist="stats/{experiment}/{sample_prefix}_map.idhist",
    params:
        sam=config["sam"],
        kmers=config["kmers"],
        mem=config["mem"],
        compression_flags=bbtools_compression_flags,
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbmap_map.benchmark.txt"
    log:
        "logs/{experiment}/bbmap/{sample_prefix}.bbmap_map.log",
    threads: 16
    shell:
        "[ -f {input.index_dir}/ref/genome/1/summary.txt ] "
        "|| (echo 'BBMap index missing summary file: {input.index_dir}/ref/genome/1/summary.txt' >&2; exit 1); "
        "[ -f {input.index_dir}/ref/genome/1/chr1.chrom.gz ] "
        "|| (echo 'BBMap index missing chrom file: {input.index_dir}/ref/genome/1/chr1.chrom.gz' >&2; exit 1); "
        "( bbmap.sh "
        "-Xmx{params.mem}g "
        "in1={input.R1_ec} "
        "in2={input.R2_ec} "
        "sam={params.sam} 32bit=t "
        "path={input.index_dir} "
        "build=1 "
        "outm=stdout.sam "
        "{params.compression_flags}"
        "k={params.kmers} "
        "t={threads} "
        "nzo=true "
        "covstats={output.covstats} "
        "basecov={output.basecov} "
        "bincov={output.bincov} "
        "ehist={output.ehist} "
        "indelhist={output.indelhist} "
        "mhist={output.mhist} "
        "idhist={output.idhist} "
        "| samtools view -b -o {output.bam} - "
        ") 2> {log}"
