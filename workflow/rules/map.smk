rule map_to_reference_bbmap:
  input:
    R1_ec=expand("results/{experiment}/{{sample_prefix}}_R1.ec.clean.trim.fastq.gz", experiment=config["experiment"]),
    R2_ec=expand("results/{experiment}/{{sample_prefix}}_R2.ec.clean.trim.fastq.gz", experiment=config["experiment"])
  output:
    mapped=temp(expand("results/{experiment}/{{sample_prefix}}.mapped.sam", experiment=config["experiment"])),
    covstats=expand("stats/{{sample_prefix}}_map.covstats"),
    covhist=expand("stats/{{sample_prefix}}_map.covhist"),
    basecov=expand("stats/{{sample_prefix}}_map.basecov"),
    bincov=expand("stats/{{sample_prefix}}_map.bincov"),
    ehist=expand("stats/{{sample_prefix}}_map.ehist"),
    indelhist=expand("stats/{{sample_prefix}}_map.indelhist"),
    mhist=expand("stats/{{sample_prefix}}_map.mhist"),
    idhist=expand("stats/{{sample_prefix}}_map.idhist")
  params:
    ref=config["reference"],
    ref_dir=config["ref_dir"],
    sam=config["sam"],
    kmers=config["kmers"],
    mem=config["mem"]
  benchmark:
    "benchmarks/{sample_prefix}.bbmap_map.benchmark.txt"
  log:
    "logs/bbmap/{sample_prefix}.bbmap_map.log"
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
    "idhist={output.idhist}"

rule sort_index_samtools:
  input:
    expand("results/{experiment}/{{sample_prefix}}.mapped.sam", experiment=config["experiment"])
  output:
    bam=temp(expand("results/{experiment}/{{sample_prefix}}.mapped.bam", experiment=config["experiment"])),
    bai=temp(expand("results/{experiment}/{{sample_prefix}}.mapped.bam.bai", experiment=config["experiment"]))
  benchmark:
    "benchmarks/{sample_prefix}.samtools_sort.benchmark.txt"
  log:
    "logs/samtools/{sample_prefix}.samtools.log"
  shell:
    "samtools sort {input} "
    " -o {output.bam} && "
    "samtools index {output.bam}"