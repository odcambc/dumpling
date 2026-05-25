rule run_rosace:
    """Run Rosace scoring for a single non-baseline condition. Snakemake
    fans this out over `experimental_conditions`; one R process per
    condition lets MCMC sampling overlap on cluster or multi-core local."""
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/enrich_format/{experiments}.tsv",
            experiments=experiment_samples,
        ),
        "results/{experiment_name}/rosace/rosace_installed.txt",
    output:
        "results/{experiment_name}/rosace/{condition}_scores.csv",
    benchmark:
        "benchmarks/{experiment_name}/{condition}.rosace.benchmark.txt"
    log:
        "logs/{experiment_name}/rosace/{condition}.rosace.log",
    threads: 4
    resources:
        mem_mb=config["mem_rosace"],
    conda:
        "../envs/rosace.yaml" if not config["rosace_local"] else None
    script:
        "scripts/run_rosace.R"


rule install_rosace:
    output:
        "results/{experiment_name}/rosace/rosace_installed.txt",
    benchmark:
        "benchmarks/{experiment_name}/install_rosace.benchmark.txt"
    log:
        "logs/{experiment_name}/install_rosace.log",
    conda:
        "../envs/rosace.yaml" if not config["rosace_local"] else None
    script:
        "scripts/install_rosace.R"
