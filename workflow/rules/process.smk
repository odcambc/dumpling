rule generate_variants:
    """Runs script to generate a designed variants csv from a DIMPLE oligo csv input."""
    input:
        oligo_file=oligo_file,
        ref_fasta=normalized_reference_file,
    output:
        expand("{variants_file}", variants_file=config["variants_file"])
        if config["regenerate_variants"]
        else [],
    benchmark:
        "benchmarks/generate_variants.benchmark.txt"
    log:
        "logs/scripts/generate_variants.log",
    script:
        "scripts/generate_variants.py"


rule process_sample:
    """Process one sample's GATK variantCounts into the per-sample Enrich2 tsv,
    processed CSV, and total-processing stats. Runs once per sample (parallel-
    safe; sibling invocations write to disjoint output paths)."""
    input:
        *([config["variants_file"]] if not config["noprocess"] else []),
        gatk="results/{experiment}/gatk/{sample_prefix}.variantCounts",
        # The designed variants file is only consumed when noprocess=False
        # (i.e. when we're filtering observed variants against the designed
        # library). Under noprocess=True it's neither read nor required to
        # exist, so don't declare it as an input then.
    output:
        enrich="results/{experiment}/processed_counts/enrich_format/{sample_prefix}.tsv",
        csv="results/{experiment}/processed_counts/{sample_prefix}.csv",
        total_stats="stats/{experiment}/processing/{sample_prefix}_total_processing.tsv",
    params:
        regenerate_variants=config["regenerate_variants"],
        gatk_dir="results/{experiment}/gatk/",
    benchmark:
        "benchmarks/{experiment}/{sample_prefix}.process_sample.benchmark.txt"
    log:
        "logs/{experiment}/scripts/{sample_prefix}.process_sample.log",
    script:
        "scripts/process_counts.py"


rule remove_zeros:
    """Runs script to remove unobserved variants from count files. This is only run
    if the remove_zeros parameter is set to True in the config file. This is used to
    handle some undesired behavior in Enrich2, and is not necessary for Rosace."""
    input:
        expand(
            "results/{{experiment}}/processed_counts/enrich_format/{samples}.tsv",
            samples=samples,
        ),
        # variants_file is not read by remove_zeros.py; it was declared as an
        # input by mistake. Dropping the declaration so the rule doesn't require
        # the designed-variants CSV when noprocess=True.
    output:
        expand(
            "results/{{experiment}}/processed_counts/removed_zeros/{samples}.tsv",
            samples=samples,
        ),
    params:
        samples=samples,
    benchmark:
        "benchmarks/{experiment}/{experiment}.remove_zeros.benchmark.txt"
    log:
        "logs/{experiment}/scripts/{experiment}.remove_zeros.log",
    script:
        "scripts/remove_zeros.py"
