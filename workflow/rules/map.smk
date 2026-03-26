rule map_to_reference_bbmap:
    """Map reads to reference sequence using BBMap. covhist output is currently disabled as it causes MultiQC bloat."""
    input:
        index_dir=f"ref/bbmap/{experiment}",
        R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
        R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
    output:
        mapped=temp("results/{experiment}/{sample_prefix}.mapped.sam"),
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
        "bbmap.sh "
        "-Xmx{params.mem}g "
        "in1={input.R1_ec} "
        "in2={input.R2_ec} "
        "sam={params.sam} 32bit=t "
        "path={input.index_dir} "
        "build=1 "
        "outm={output.mapped} "
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
        "idhist={output.idhist} 2> {log}"


rule sort_index_samtools:
    """Sort and index mapped reads using samtools."""
    input:
        "results/{experiment}/{sample_prefix}.mapped.sam",
    output:
        bam=temp("results/{experiment}/{sample_prefix}.mapped.bam"),
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.samtools_sort.benchmark.txt"
    log:
        "logs/{experiment}/samtools/{sample_prefix}.samtools.log",
    params:
        extra="-O bam --write-index -o {output.bam}",
    wrapper:
        "v3.1.0/bio/samtools/sort"
