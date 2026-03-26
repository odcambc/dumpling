rule gatk_ASM:
    """Use the GATK AnalyzeSaturationMutagenesis tool to call and count variants.
    This is the culmination of the initial variant calling pipeline."""
    input:
        index=str(reference_file) + ".fai",
        ref_dict=str(reference_file.parent / f"{reference_name}.dict"),
        sam="results/{experiment}/{sample_prefix}.mapped.sam",
    output:
        "results/{experiment}/gatk/{sample_prefix}.variantCounts",
    params:
        orf=config["orf"],
        min_q=config["min_q"],
        min_variant_obs=config["min_variant_obs"],
        ref=str(reference_file),
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
        "--min-variant-obs {params.min_variant_obs} "
        "--dont-ignore-disjoint-pairs "
        "-O results/{wildcards.experiment}/gatk/{wildcards.sample_prefix} 2>{log}"
