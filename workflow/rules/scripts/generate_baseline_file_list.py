import logging
from script_utils import load_experiments, run_script


"""Generate a list of files specific to baseline conditionsto be used as input for MultiQC.
This is just a list of files that will be passed using --file-list

The baseline QC reads the following files:
- fastqc
- bbduk
- bbmap
- gatk analyzesaturationmutagenesis
- count processing script

Essentially, it does not continue to Enrich2, since we don't run that on baseline samples.
"""


def _run(snakemake):
    output_file = snakemake.output[0]
    experiment_name = snakemake.config["experiment"]
    fastqc_names = snakemake.params["fastqc_names"]
    # The map-stage stats set depends on the aligner: bbmap emits the
    # BBTools-format histograms; minimap2 emits samtools-format outputs.
    # The map_to_reference_minimap2 rule (workflow/rules/map.smk) and the
    # multiqc_dir rule (workflow/rules/qc.smk) follow the same convention.
    aligner = snakemake.config.get("aligner", "bbmap")

    logging.debug("Writing to file: %s", output_file)

    baseline_condition = snakemake.config["baseline_condition"]

    experiments = load_experiments(snakemake.config["experiment_file"])

    baseline_samples = experiments.loc[
        experiments["condition"] == baseline_condition, "sample"
    ]

    stats_prefix = f"./stats/{experiment_name}"
    results_prefix = f"./results/{experiment_name}"

    # Upstream trim/clean/correct + bbmerge stats are aligner-independent.
    upstream_stats_exts = (
        "_trim.qhist",
        "_trim.bhist",
        "_trim.gchist",
        "_trim.aqhist",
        "_trim.lhist",
        "_trim.stats.txt",
        "_trim_contam.stats.txt",
        "_merge.ihist",
    )

    if aligner == "bbmap":
        map_stats_exts = (
            "_map.covstats",
            "_map.basecov",
            "_map.bincov",
            "_map.ehist",
            "_map.indelhist",
            "_map.mhist",
            "_map.idhist",
        )
    else:  # minimap2
        # samtools coverage was intentionally dropped: it does per-base
        # pileup and costs ~14s/sample on deep DMS data. Re-add via
        # mosdepth in a follow-up — tracked in tasks/tasks.md.
        map_stats_exts = (
            "_samtools_stats.txt",
            "_samtools_flagstat.txt",
        )

    with open(output_file, "w+") as f:
        logging.debug("Writing to file: %s", output_file)

        for sample in baseline_samples.index:
            # FastQC output is named from the *fastq filename*, not the
            # sample name from the experiment CSV. Look up the resolved
            # base names via the same fastqc_names map the producing rule
            # uses (computed in common.smk).
            file_prefix = experiments.loc[sample, "file"]
            fastqc_R1 = fastqc_names[file_prefix]["R1"]
            fastqc_R2 = fastqc_names[file_prefix]["R2"]
            f.write(f"{stats_prefix}/fastqc/{fastqc_R1}_fastqc.html\n")
            f.write(f"{stats_prefix}/fastqc/{fastqc_R2}_fastqc.html\n")

            # BBDuk / BBMerge / GATK rules emit per-sample-name outputs
            # (rule wildcards are `{sample_prefix}` = the sample name), so
            # these path templates are correct as-is.
            for ext in upstream_stats_exts + map_stats_exts:
                f.write(f"{stats_prefix}/{sample}{ext}\n")

            # GATK ASM side-output stats files (live next to the .variantCounts).
            for ext in (".coverageLengthCounts", ".readCounts", ".refCoverage"):
                f.write(f"{results_prefix}/gatk/{sample}{ext}\n")

            # Processed counts CSV. process_counts.py writes to
            # `processed_counts/{sample}.csv` — no `counts/` subdirectory.
            f.write(f"{results_prefix}/processed_counts/{sample}.csv\n")

    logging.debug("Done")


def main():
    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
