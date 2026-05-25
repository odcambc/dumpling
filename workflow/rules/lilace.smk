rule run_lilace:
    """Run Lilace scoring for a single non-baseline condition. Snakemake
    fans this out over `experimental_conditions`; one R process per
    condition lets Stan sampling overlap on cluster or multi-core local."""
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/enrich_format/{experiments}.tsv",
            experiments=experiment_samples,
        ),
        "results/{experiment_name}/lilace/lilace_installed.txt",
    output:
        "results/{experiment_name}/lilace/{condition}_scores.csv",
    benchmark:
        "benchmarks/{experiment_name}/{condition}.lilace.benchmark.txt"
    log:
        "logs/{experiment_name}/lilace/{condition}.lilace.log",
    threads: 4
    resources:
        mem_mb=config["mem_lilace"],
    conda:
        "../envs/lilace.yaml" if not config["lilace_local"] else None
    script:
        "scripts/run_lilace.R"


rule install_lilace:
    output:
        "results/{experiment_name}/lilace/lilace_installed.txt",
    benchmark:
        "benchmarks/{experiment_name}/install_lilace.benchmark.txt"
    log:
        "logs/{experiment_name}/install_lilace.log",
    conda:
        "../envs/lilace.yaml" if not config["lilace_local"] else None
    script:
        "scripts/install_lilace.R"
