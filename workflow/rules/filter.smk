rule bbduk_trim_adapters:
  input:
    R1=expand("{data_dir}/{{sample_prefix}}_R1_001.fastq.gz", data_dir=config["data_dir"]),
    R2=expand("{data_dir}/{{sample_prefix}}_R2_001.fastq.gz", data_dir=config["data_dir"])
  output:
    R1_trim=temp("results/{experiment}/{sample_prefix}_R1.trim.fastq.gz"),
    R2_trim=temp("results/{experiment}/{sample_prefix}_R2.trim.fastq.gz"),
    qhist="stats/{experiment}/{sample_prefix}_trim.qhist",
    bhist="stats/{experiment}/{sample_prefix}_trim.bhist",
    gchist="stats/{experiment}/{sample_prefix}_trim.gchist",
    aqhist="stats/{experiment}/{sample_prefix}_trim.aqhist",
    lhist="stats/{experiment}/{sample_prefix}_trim.lhist",
    stats="stats/{experiment}/{sample_prefix}_trim.stats.txt"
  params:
    adapters=config["adapters"],
    kmers=config["kmers"]
  benchmark:
    "benchmarks/{experiment}/{sample_prefix}.bbduk_trim.benchmark.txt"
  log:
    "logs/{experiment}/bbduk/{sample_prefix}.trim.bbduk.log"
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
    "gcbins=auto 2> {log}"

rule remove_contaminants:
  input:
    R1_trim="results/{experiment}/{sample_prefix}_R1.trim.fastq.gz",
    R2_trim="results/{experiment}/{sample_prefix}_R2.trim.fastq.gz"
  output:
    R1_clean=temp("results/{experiment}/{sample_prefix}_R1.clean.trim.fastq.gz"),
    R2_clean=temp("results/{experiment}/{sample_prefix}_R2.clean.trim.fastq.gz"),
    stats="stats/{experiment}/{sample_prefix}_trim.stats.txt"
  params:
    contaminant_ref=config["contaminant_ref"]
  benchmark:
    "benchmarks/{experiment}/{sample_prefix}.bbduk_clean.benchmark.txt"
  log:
    "logs/{experiment}/bbduk/{sample_prefix}.clean.bbduk.log"
  shell:
    "bbduk.sh "
    "in={input.R1_trim} "
    "in2={input.R2_trim} "
    "out={output.R1_clean} "
    "out2={output.R2_clean} "
    "stats={output.stats} "
    "overwrite=true "
    "k=31 "
    "ref={params.contaminant_ref} 2> {log}"

rule error_correct_bbmerge:
  input:
    R1_clean="results/{experiment}/{sample_prefix}_R1.clean.trim.fastq.gz",
    R2_clean="results/{experiment}/{sample_prefix}_R2.clean.trim.fastq.gz"
  output:
    R1_ec=temp("results/{experiment}/{sample_prefix}_R1.ec.clean.trim.fastq.gz"),
    R2_ec=temp("results/{experiment}/{sample_prefix}_R2.ec.clean.trim.fastq.gz"),
    ihist="stats/{experiment}/{sample_prefix}_merge.ihist"
  params:
  benchmark:
    "benchmarks/{experiment}/{sample_prefix}.bbmerge.benchmark.txt"
  log:
    "logs/{experiment}/bbmerge/{sample_prefix}.bbmerge.log"
  shell:
    "bbmerge.sh "
    "in={input.R1_clean} "
    "in2={input.R2_clean} "
    "out={output.R1_ec} "
    "out2={output.R2_ec} "
    "overwrite=true "
    "ihist={output.ihist} "
    "ecco mix showhiststats=t 2> {log}"