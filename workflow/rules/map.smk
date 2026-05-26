if config["aligner"] == "bbmap":

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

elif config["aligner"] == "minimap2":

    rule map_to_reference_minimap2:
        """Map reads with minimap2 in short-read mode (`-ax sr`), streaming SAM
        through `samtools view -b` so the on-disk intermediate is a BAM. Then
        run `samtools stats` and `samtools flagstat` on that BAM as the
        per-sample mapping-stage QC artifacts.

        minimap2 doesn't emit the BBTools-format histograms BBMap does
        (`_map.covstats`, `_map.ehist`, etc.), but MultiQC's samtools-stats /
        samtools-flagstat parsers cover most load-bearing DMS QC: mapping
        rate, insert sizes, global error rate, indel size distribution,
        read length, GC content. Mean depth can be derived from
        `bases mapped (cigar)` ÷ `reference length` in the stats output.

        Per-position coverage (`samtools coverage` / `samtools depth`) is
        explicitly NOT computed here: both do per-base pileup and cost
        ~14s/sample on the example fixture (deeply-covered DMS data is
        ~60,000× per position, so the per-base-event work dominates).
        Real DMS prod depths are similar or higher, so this isn't a fixture
        artifact. **TODO**: re-add coverage via mosdepth (purpose-built
        for fast deep-coverage; ~10× faster than samtools coverage; has
        dedicated MultiQC parser). Tracked in `tasks/tasks.md`.

        Downstream rules that depend on map-stage stats branch on
        `config['aligner']`: `multiqc_dir` (qc.smk) and the baseline file
        list generator (`generate_baseline_file_list.py`) both pick up the
        samtools paths under `aligner: minimap2`."""
        input:
            index=f"ref/minimap2/{experiment}_{ref_digest}.mmi",
            R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
            R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
        output:
            bam=temp("results/{experiment}/{sample_prefix}.mapped.bam"),
            stats="stats/{experiment}/{sample_prefix}_samtools_stats.txt",
            flagstat="stats/{experiment}/{sample_prefix}_samtools_flagstat.txt",
        benchmark:
            "benchmarks/{experiment}/{sample_prefix}.minimap2_map.benchmark.txt"
        log:
            "logs/{experiment}/minimap2/{sample_prefix}.minimap2_map.log",
        threads: 16
        shell:
            # samtools stats and flagstat both scan the BAM linearly without
            # caring about sort order, so the BAM stays unsorted (matches
            # the bbmap rule's behaviour and what GATK ASM expects).
            "( minimap2 -ax sr --MD -t {threads} "
            "{input.index} {input.R1_ec} {input.R2_ec} "
            "| samtools view -b -o {output.bam} - "
            ") 2> {log}; "
            "samtools stats    -@ {threads} {output.bam} > {output.stats}; "
            "samtools flagstat -@ {threads} {output.bam} > {output.flagstat}"
