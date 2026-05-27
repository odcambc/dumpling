rule format_mavedb:
    """Format variant scores into MaveDB score set CSV format.

    Run for a specific condition:
      snakemake results/<experiment>/deposit/mavedb/<condition>_mavedb.csv
    """
    input:
        scores=f"results/{{experiment_name}}/{scoring_backend}/{{condition}}_scores.csv",
    output:
        "results/{experiment_name}/deposit/mavedb/{condition}_mavedb.csv",
    params:
        backend=scoring_backend,
    log:
        "logs/{experiment_name}/deposit/{condition}_mavedb.log",
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
    params:
        data_dir=config["data_dir"],
        sra_config=config.get("sra", {}),
    log:
        "logs/{experiment_name}/deposit/sra.log",
    script:
        "scripts/prepare_sra.py"
