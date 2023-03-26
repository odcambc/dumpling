rule map_to_reference_bbmap:
  input:
    R1_ec="results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz",
    R2_ec="results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz"
  output:
    mapped=temp("results/{experiment}/{sample_prefix}.mapped.sam"),
    covstats="stats/{experiment}/{sample_prefix}_map.covstats",
    covhist="stats/{experiment}/{sample_prefix}_map.covhist",
    basecov="stats/{experiment}/{sample_prefix}_map.basecov",
    bincov="stats/{experiment}/{sample_prefix}_map.bincov",
    ehist="stats/{experiment}/{sample_prefix}_map.ehist",
    indelhist="stats/{experiment}/{sample_prefix}_map.indelhist",
    mhist="stats/{experiment}/{sample_prefix}_map.mhist",
    idhist="stats/{experiment}/{sample_prefix}_map.idhist"
  params:
    ref=config["reference"],
    ref_dir=config["ref_dir"],
    sam=config["sam"],
    kmers=config["kmers"],
    mem=config["mem"]
  benchmark:
    "benchmarks/{experiment}/{sample_prefix}.bbmap_map.benchmark.txt"
  log:
    "logs/{experiment}/bbmap/{sample_prefix}.bbmap_map.log"
  threads: 
    16
  shell:
    "bbmap.sh "
    "-Xmx{params.mem}g "
    "in1={input.R1_ec} "
    "in2={input.R2_ec} "
    "sam={params.sam} 32bit=t "
    "ref={params.ref_dir}/{params.ref} "
    "outm={output.mapped} "
    "k={params.kmers} "
    "t={threads} "
    "covstats={output.covstats} "
    "covhist={output.covhist} "
    "basecov={output.basecov} "
    "bincov={output.bincov} "
    "ehist={output.ehist} "
    "indelhist={output.indelhist} "
    "mhist={output.mhist} "
    "idhist={output.idhist} 2> {log}"

rule sort_index_samtools:
  input:
    "results/{experiment}/{sample_prefix}.mapped.sam"
  output:
    bam=temp("results/{experiment}/{sample_prefix}.mapped.bam"),
    bai=temp("results/{experiment}/{sample_prefix}.mapped.bam.bai")
  benchmark:
    "benchmarks/{experiment}/{sample_prefix}.samtools_sort.benchmark.txt"
  log:
    "logs/{experiment}/samtools/{sample_prefix}.samtools.log"
  shell:
    "samtools sort {input} "
    " -o {output.bam} && "
    "samtools index {output.bam} 2> {log}"