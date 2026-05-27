def _mavedb_condition_samples(condition):
    """Return the ordered list of sample names whose processed_counts/{sample}.csv
    should appear as columns in {condition}_mavedb_counts.csv.

    This explicitly only returns conditions referring to _experimental conditions_
    since those are the basis for any scores. The expectation is that full
    library QC is deposited elsewhere, given baselines are optional.

    Not that by construction of the dag and workflow "condition" is solely
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
