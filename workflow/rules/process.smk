rule process_counts:
    """Runs script to process variant counts generated by GATK ASM into a table of counts per sample, 
    as is required for running Enrich2 and Rosace."""
    input:
        expand("results/{{experiment}}/gatk/{samples}.variantCounts", samples=samples),
    output:
        expand("results/{{experiment}}/processed_counts/{samples}.tsv", samples=samples),
        expand(
            "stats/{{experiment}}/processing/{samples}_total_processing.tsv",
            samples=samples,
        ),
    params:
        remove_zeros=config["remove_zeros"],
        regenerate_variants=config["regenerate_variants"],
        samples=samples,
        gatk_dir="results/{experiment}/gatk/",
    benchmark:
        "benchmarks/{experiment}/{experiment}.process_counts.benchmark.txt"
    log:
        "logs/{experiment}/scripts/{experiment}.process_counts.log",
    script:
        "scripts/process_counts.py"
