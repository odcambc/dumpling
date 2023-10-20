rule bbduk_trim_adapters:
    """Remove adapters with BBDuk. This uses all of the adapters in the adapters.fa file, which is included with BBTools.
    This also performs the first sample renaming step, using the mapping file -> sample from the provided experiment CSV."""
    input:
        unpack(get_file_from_sample),
    output:
        R1_trim="results/{experiment}/{sample_prefix}_R1.trim.fastq.gz",
        R2_trim="results/{experiment}/{sample_prefix}_R2.trim.fastq.gz",
        qhist="stats/{experiment}/{sample_prefix}_trim.qhist",
        bhist="stats/{experiment}/{sample_prefix}_trim.bhist",
        gchist="stats/{experiment}/{sample_prefix}_trim.gchist",
        aqhist="stats/{experiment}/{sample_prefix}_trim.aqhist",
        lhist="stats/{experiment}/{sample_prefix}_trim.lhist",
        stats="stats/{experiment}/{sample_prefix}_trim.stats.txt",
    params:
        adapters=adapters_ref
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbduk_trim.benchmark.txt"
    log:
        "logs/{experiment}/bbduk/{sample_prefix}.trim.bbduk.log",
    threads: 16
    shell:
        "bbduk.sh in1={input.R1} in2={input.R2} "
        "ref={params.adapters} ktrim=r k=23 mink=10 hdist=1 ftr=149 tpe tbo "
        "out1={output.R1_trim} "
        "out2={output.R2_trim} "
        "bhist={output.bhist} "
        "qhist={output.qhist} "
        "gchist={output.gchist} "
        "aqhist={output.aqhist} "
        "lhist={output.lhist} "
        "stats={output.stats} "
        "overwrite=true "
        "t={threads} "
        "gcbins=auto 2> {log}"


rule remove_contaminants:
    """Remove typical contaminants with BBDuk.
    The contaminants files used, by default, are included with BBTools, and consist of PhiX and some (possibly JGI-specific) other sequencing artifacts."""
    input:
        R1_trim="results/{experiment}/{sample_prefix}_R1.trim.fastq.gz",
        R2_trim="results/{experiment}/{sample_prefix}_R2.trim.fastq.gz",
    output:
        R1_clean="results/{experiment}/{sample_prefix}_R1.clean.trim.fastq.gz",
        R2_clean="results/{experiment}/{sample_prefix}_R2.clean.trim.fastq.gz",
        stats="stats/{experiment}/{sample_prefix}_trim.stats.txt",
    params:
        prefix=os.environ["CONDA_PREFIX"],
        contaminants=contaminants_ref,
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbduk_clean.benchmark.txt"
    log:
        "logs/{experiment}/bbduk/{sample_prefix}.clean.bbduk.log",
    threads: 16
    shell:
        "bbduk.sh "
        "in={input.R1_trim} "
        "in2={input.R2_trim} "
        "out={output.R1_clean} "
        "out2={output.R2_clean} "
        "stats={output.stats} "
        "overwrite=true "
        "k=31 "
        "t={threads} "
        "ref={params.contaminants} 2> {log}"


rule error_correct_bbmerge:
    """Perform error correction using paired-end read overlap using BBMerge.
    This keeps separate R1 and R2 files and does not actually merge the reads."""
    input:
        R1_clean="results/{experiment}/{sample_prefix}_R1.clean.trim.fastq.gz",
        R2_clean="results/{experiment}/{sample_prefix}_R2.clean.trim.fastq.gz",
    output:
        R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
        R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz",
        ihist="stats/{experiment}/{sample_prefix}_merge.ihist",
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.bbmerge.benchmark.txt"
    log:
        "logs/{experiment}/bbmerge/{sample_prefix}.bbmerge.log",
    threads: 16
    shell:
        "bbmerge.sh "
        "in={input.R1_clean} "
        "in2={input.R2_clean} "
        "out={output.R1_ec} "
        "out2={output.R2_ec} "
        "overwrite=true "
        "ihist={output.ihist} "
        "t={threads} "
        "ecco mix showhiststats=t 2> {log}"
