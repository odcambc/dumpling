rule bbduk_trim_adapters:
  input:
    R1=expand("{data_dir}/{experiment_dir}/{{sample_prefix}}_R1_001.fastq.gz", data_dir=config["data_dir"], experiment_dir=config["experiment_dir"]),
    R2=expand("{data_dir}/{experiment_dir}/{{sample_prefix}}_R2_001.fastq.gz", data_dir=config["data_dir"], experiment_dir=config["experiment_dir"])
  output:
    R1_trim=temp(expand("results/{experiment}/{{sample_prefix}}_R1.trim.fastq.gz", experiment=config["experiment"])),
    R2_trim=temp(expand("results/{experiment}/{{sample_prefix}}_R2.trim.fastq.gz", experiment=config["experiment"])),
    qhist=expand("stats/{{sample_prefix}}_trim.qhist"),
    bhist=expand("stats/{{sample_prefix}}_trim.bhist"),
    gchist=expand("stats/{{sample_prefix}}_trim.gchist"),
    aqhist=expand("stats/{{sample_prefix}}_trim.aqhist"),
    lhist=expand("stats/{{sample_prefix}}_trim.lhist"),
    stats=expand("stats/{{sample_prefix}}_trim.stats.txt")
  params:
    adapters=config["adapters"],
    kmers=config["kmers"]
  benchmark:
      "benchmarks/{sample_prefix}.bbduk_trim.benchmark.txt"
  log:
    "logs/bbduk/{sample_prefix}.trim.bbduk.log"
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
    R1_trim=expand("results/{experiment}/{{sample_prefix}}_R1.trim.fastq.gz", experiment=config["experiment"]),
    R2_trim=expand("results/{experiment}/{{sample_prefix}}_R2.trim.fastq.gz", experiment=config["experiment"])
  output:
    R1_clean=temp(expand("results/{experiment}/{{sample_prefix}}_R1.clean.trim.fastq.gz", experiment=config["experiment"])),
    R2_clean=temp(expand("results/{experiment}/{{sample_prefix}}_R2.clean.trim.fastq.gz", experiment=config["experiment"])),
    stats=expand("stats/{{sample_prefix}}_trim.stats.txt")
  params:
    contaminant_ref=config["contaminant_ref"]
  benchmark:
    "benchmarks/{sample_prefix}.bbduk_clean.benchmark.txt"
  log:
    "logs/bbduk/{sample_prefix}.clean.bbduk.log"
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
    R1_clean=expand("results/{experiment}/{{sample_prefix}}_R1.clean.trim.fastq.gz", experiment=config["experiment"]),
    R2_clean=expand("results/{experiment}/{{sample_prefix}}_R2.clean.trim.fastq.gz", experiment=config["experiment"])
  output:
    R1_ec=temp(expand("results/{experiment}/{{sample_prefix}}_R1.ec.clean.trim.fastq.gz", experiment=config["experiment"])),
    R2_ec=temp(expand("results/{experiment}/{{sample_prefix}}_R2.ec.clean.trim.fastq.gz", experiment=config["experiment"])),
    ihist="stats/{sample_prefix}_merge.ihist"
  params:
  benchmark:
    "benchmarks/{sample_prefix}.bbmerge.benchmark.txt"
  log:
    "logs/bbmerge/{sample_prefix}.bbmerge.log"
  shell:
    "bbmerge.sh "
    "in={input.R1_clean} "
    "in2={input.R2_clean} "
    "out={output.R1_ec} "
    "out2={output.R2_ec} "
    "overwrite=true "
    "ihist={output.ihist} "
    "ecco mix showhiststats=t 2>{log}"