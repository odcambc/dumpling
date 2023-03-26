rule gatk_ASM:
  input:
    index=expand("{ref_dir}/{ref}.fasta.fai", ref_dir=config["ref_dir"], ref=config["reference2"]),
    ref_dict=expand("{ref_dir}/{ref}.dict", ref_dir=config["ref_dir"], ref=config["reference2"]),
    sam=expand("results/{experiment}/{{sample_prefix}}.mapped.sam", experiment=config["experiment"])
  output:
    expand("results/{experiment}/gatk/{{sample_prefix}}.aaCounts", experiment=config["experiment"])
  params:
    orf=config["orf"],
    kmers=config["kmers"],
    ref=expand("{ref_dir}/{ref}.fasta", ref_dir=config["ref_dir"], ref=config["reference2"]),
  benchmark: "benchmarks/{sample_prefix}.gatk_asm.benchmark.txt"
  log: "logs/gatk/{sample_prefix}.gatk_asm.log"
  default_target: True
  threads: 1
  shell:
    "gatk AnalyzeSaturationMutagenesis "
    "-I {input.sam} "
    "-R {params.ref} "
    "--orf {params.orf} "
    "-O results/{config[experiment]}/gatk/{wildcards.sample_prefix} 2>{log}"
    