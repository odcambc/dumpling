rule run_lilace:
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/enrich_format/{experiments}.tsv",
            experiments=experiment_samples,
        ),
        "results/{experiment_name}/lilace/lilace_installed.txt",
    output:
        expand(
            "results/{{experiment_name}}/lilace/{conditions}_scores.csv",
            conditions=experimental_conditions,
        ),
    benchmark:
        "benchmarks/{experiment_name}/{experiment_name}.lilace.benchmark.txt"
    log:
        "logs/{experiment_name}/lilace/{experiment_name}.lilace.log",
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
