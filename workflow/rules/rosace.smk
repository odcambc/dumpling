rule run_rosace:
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/{experiments}.tsv",
            experiments=experiment_samples,
        ),
    output:
        expand(
            "results/{{experiment_name}}/rosace/scores_{conditions}.tsv",
            conditions=experimental_conditions,
        ),
    benchmark:
        "benchmarks/{experiment_name}/{experiment_name}.rosace.benchmark.txt"
    log:
        "logs/{experiment_name}/rosace/{experiment_name}.rosace.log",
    script:
        "scripts/run_rosace.R"
