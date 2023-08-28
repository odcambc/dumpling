rule process_counts:
    input:
        expand(
            "results/{{experiment}}/gatk/{sample_name}.variantCounts",
            condition_samples=config["condition_samples"],
        ),
    output:
        "results/counts/{experiment}/{sample_name}.csv",
    params:
        remove_zeros=config["remove_zeros"],
        regenerate_variants=config["regenerate_variants"],
    benchmark:
        "benchmarks/{experiment}/{sample_name}.process_counts.benchmark.txt"
    log:
        "logs/{experiment}/process_counts/{sample_name}.trim.bbduk.log",
    script:
        "scripts/process_counts.py"
    shell:
        "python {script} remove_zeros"


#    expand("results/{experiment_name}/gatk/{sample_name}.variantCounts", sample_name=config["samples"], experiment_name=config["experiment"]),
