def _mavedb_condition_samples(condition):
    """Return the ordered list of sample names whose processed_counts/{sample}.csv
    should appear as columns in {condition}_mavedb_counts.csv.

    This explicitly only returns conditions referring to _experimental conditions_
    since those are the basis for any scores. The expectation is that full
    library QC is deposited elsewhere, given baselines are optional.

    Note that by construction of the dag and workflow "condition" is solely
    experimental conditions and so we simply read the row corresponding to that.
    """
    experimental_condition_rows = experiments[
        experiments["condition"] == condition
    ].sort_values(["replicate", "time"])

    return experimental_condition_rows["sample"].tolist()


def _mavedb_condition_count_files(wildcards):
    samples = _mavedb_condition_samples(wildcards.condition)
    return expand(
        "results/{experiment_name}/processed_counts/{sample}.csv",
        experiment_name=wildcards.experiment_name,
        sample=samples,
    )


rule format_mavedb:
    """Format variant scores and raw counts into MaveDB score/count CSVs.

Run for a specific condition:
snakemake results/<experiment>/deposit/mavedb/<condition>_mavedb.csv
snakemake results/<experiment>/deposit/mavedb/<condition>_mavedb_counts.csv

Both outputs share the same variant set (union of scored + counted
variants) and the same hgvs_pro index column, per the MaveDB spec
requirement that score and count tables align row-for-row.
"""
    input:
        scores=f"results/{{experiment_name}}/{scoring_backend}/{{condition}}_scores.csv",
        counts=_mavedb_condition_count_files,
    output:
        scores="results/{experiment_name}/deposit/mavedb/{condition}_mavedb.csv",
        counts="results/{experiment_name}/deposit/mavedb/{condition}_mavedb_counts.csv",
    log:
        "logs/{experiment_name}/deposit/{condition}_mavedb.log",
    params:
        backend=scoring_backend,
        sample_names=lambda wc: _mavedb_condition_samples(wc.condition),
    script:
        "scripts/format_mavedb.py"


rule format_cosmos:
    """Format per-condition variant scores into a single cosmos input CSV.

Run with:
snakemake results/<experiment>/cosmos/<experiment>_cosmos.csv

Unlike format_mavedb (one file per condition), this is multi-condition: it
joins the score CSVs of every condition assigned a `phenotype` slot in the
experiment CSV into one wide table with beta_hat_N/se_hat_N column pairs, in
slot order. See docs/cosmos_export_design.md.
"""
    input:
        # Score CSVs in slot order: cosmos_phenotype_conditions[i] -> beta_hat_{i+1}.
        scores=expand(
            f"results/{{experiment_name}}/{scoring_backend}/{{condition}}_scores.csv",
            condition=cosmos_phenotype_conditions,
            allow_missing=True,
        ),
    output:
        cosmos="results/{experiment_name}/cosmos/{experiment_name}_cosmos.csv",
    log:
        "logs/{experiment_name}/cosmos.log",
    params:
        backend=scoring_backend,
    script:
        "scripts/format_cosmos.py"


rule run_cosmos:
    """Run cosmos on the formatted input: the per-position direct/indirect
effect decomposition across the two phenotype conditions (slots 1 -> 2).

Run with:
snakemake results/<experiment>/cosmos/<experiment>_cosmos_results.csv

NOTE: cosmos fits a model PER POSITION (~tens of seconds each) and this runs
them serially, so a large library takes hours. Parallelizing across positions
is a tracked future optimization (tasks.md). See docs/cosmos_export_design.md.
"""
    input:
        cosmos="results/{experiment_name}/cosmos/{experiment_name}_cosmos.csv",
    output:
        results="results/{experiment_name}/cosmos/{experiment_name}_cosmos_results.csv",
    log:
        "logs/{experiment_name}/cosmos_run.log",
    params:
        # The two phenotype conditions in slot order (slot 1 -> beta_hat_1 ->
        # cosmos x; slot 2 -> beta_hat_2 -> cosmos y). Labels only.
        phenotypes=cosmos_phenotype_conditions,
    resources:
        mem_mb=config["mem_cosmos"],
    conda:
        "../envs/cosmos.yaml"
    script:
        "scripts/run_cosmos.py"


rule prepare_sra:
    """Prepare SRA submission metadata and FASTQ file list.

Run with:
snakemake results/<experiment>/deposit/sra/sra_metadata.tsv
"""
    input:
        experiment_file=config["experiment_file"],
    output:
        metadata="results/{experiment_name}/deposit/sra/sra_metadata.tsv",
        filelist="results/{experiment_name}/deposit/sra/sra_files.txt",
    log:
        "logs/{experiment_name}/deposit/sra.log",
    params:
        data_dir=config["data_dir"],
        sra_config=config.get("sra", {}),
    script:
        "scripts/prepare_sra.py"
