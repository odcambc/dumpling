rule gatk_ASM:
  input:
    index=expand("{ref_dir}/{ref}.fasta.fai", ref_dir=config["ref_dir"], ref=config["reference2"]),
    ref_dict=expand("{ref_dir}/{ref}.dict", ref_dir=config["ref_dir"], ref=config["reference2"]),
    sam="results/{experiment}/{sample_prefix}.mapped.sam"
  output:
    "results/{experiment}/gatk/{sample_prefix}.aaCounts"
  params:
    orf=config["orf"],
    kmers=config["kmers"],
    ref=expand("{ref_dir}/{ref}.fasta", ref_dir=config["ref_dir"], ref=config["reference2"]),
  benchmark: "benchmarks/{experiment}/{sample_prefix}.gatk_asm.benchmark.txt"
  log: "logs/{experiment}/{sample_prefix}.gatk_asm.benchmark.txt"
  default_target: True
  threads: 1
  shell:
    "gatk AnalyzeSaturationMutagenesis "
    "-I {input.sam} "
    "-R {params.ref} "
    "--orf {params.orf} "
    "-O results/{wildcards.experiment}/gatk/{wildcards.sample_prefix} 2>{log}"
    