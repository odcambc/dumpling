rule gatk_ASM:
    """Use the GATK AnalyzeSaturationMutagenesis tool to call and count variants.
    This is the culmination of the initial variant calling pipeline."""
    input:
        index=expand(
            "{ref_dir}/{ref_name}.fasta.fai",
            ref_dir=config["ref_dir"],
            ref_name=reference_name,
        ),
        ref_dict=expand(
            "{ref_dir}/{ref_name}.dict",
            ref_dir=config["ref_dir"],
            ref_name=reference_name,
        ),
        sam="results/{experiment}/{sample_prefix}.mapped.sam",
    output:
        "results/{experiment}/gatk/{sample_prefix}.variantCounts",
    params:
        orf=config["orf"],
        min_q=config["min_q"],
        ref=expand(
            "{ref_dir}/{ref_name}.fasta",
            ref_dir=config["ref_dir"],
            ref_name=reference_name,
        ),
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.gatk_asm.benchmark.txt"
    log:
        "logs/{experiment}/gatk/{sample_prefix}.gatk_asm.log",
    threads: 1
    shell:
        "gatk AnalyzeSaturationMutagenesis "
        "-I {input.sam} "
        "-R {params.ref} "
        "--orf {params.orf} "
        "--sequence-dictionary {input.ref_dict} "
        "--min-q {params.min_q} "
        "-O results/{wildcards.experiment}/gatk/{wildcards.sample_prefix} 2>{log}"
