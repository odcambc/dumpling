rule run_rosace_aa:
    """Run rosace-aa scoring for a single non-baseline condition. Mirrors
    `run_rosace`'s shape: Snakemake fans this rule out over
    `experimental_conditions`, one R process per condition. rosace-aa extends
    rosace with position + amino-acid substitution effect decomposition;
    its score CSV layout matches rosace's (so format_mavedb works unchanged),
    but the entrypoint takes additional metadata-column args.

    See workflow/rules/scripts/run_rosace_aa.R for the RunRosace call site
    and the wt.col/mut.col/aa.code/stop.col/stop.name parameter mapping
    (audit notes in tasks/tasks.md)."""
    input:
        expand(
            "results/{{experiment_name}}/processed_counts/enrich_format/{experiments}.tsv",
            experiments=experiment_samples,
        ),
        "results/{experiment_name}/rosace_aa/rosace_aa_installed.txt",
    output:
        "results/{experiment_name}/rosace_aa/{condition}_scores.csv",
    benchmark:
        "benchmarks/{experiment_name}/{condition}.rosace_aa.benchmark.txt"
    log:
        "logs/{experiment_name}/rosace_aa/{condition}.rosace_aa.log",
    threads: 4
    resources:
        mem_mb=config["mem_rosace_aa"],
    conda:
        "../envs/rosace_aa.yaml" if not config["rosace_aa_local"] else None
    script:
        "scripts/run_rosace_aa.R"


rule install_rosace_aa:
    output:
        "results/{experiment_name}/rosace_aa/rosace_aa_installed.txt",
    benchmark:
        "benchmarks/{experiment_name}/install_rosace_aa.benchmark.txt"
    log:
        "logs/{experiment_name}/install_rosace_aa.log",
    conda:
        "../envs/rosace_aa.yaml" if not config["rosace_aa_local"] else None
    script:
        "scripts/install_rosace_aa.R"
