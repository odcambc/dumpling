rule run_rosace:
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/{experiments}.tsv",
            experiments=experiment_samples,
        ),
    output:
        expand(
            "results/{{experiment_name}}/rosace/{conditions}_scores.tsv",
            conditions=experimental_conditions,
        ),
    benchmark:
        "benchmarks/{experiment_name}/{experiment_name}.rosace.benchmark.txt"
    log:
        "logs/{experiment_name}/rosace/{experiment_name}.rosace.log",
    threads: 4
    script:
        "scripts/run_rosace.R"

rule install_rosace:
    output:
        "renv/activate.R"
    benchmark:
        "benchmarks/install_rosace.benchmark.txt"
    log:
        "logs/install_rosace.log",
    script:
        "scripts/install_rosace.R"
