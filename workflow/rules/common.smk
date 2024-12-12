def get_file_from_sample(wildcards):
    """Maps from the sequencing output file names to the sample names defined in the experiment CSV.
    This is used in rule bbduk_trim_adapters"""
    filename = experiments.loc[
        experiments["sample"] == wildcards.sample_prefix, "file"
    ].squeeze()
    prefix = config["data_dir"] + "/" + filename
    R1 = prefix + "_R1_001.fastq.gz"
    R2 = prefix + "_R2_001.fastq.gz"
    return {"R1": R1, "R2": R2}


def get_ref(wildcards):
    """Removes file suffix from reference fasta file. This is used in rule bwa_index"""
    prefix = config["reference"].split(".fasta")[0]
    return prefix


def pass_names(names):
    """Takes either an array of strings (such as files) or a single string and returns a comma-separated string.
    This is used while passing references to bbduk."""
    if isinstance(names, str):
        return names
    else:
        return ",".join(names)


def get_baseline_samples(experiments, samples):
    """Returns a list of baseline sample and files"""
    baseline_samples = []
    baseline_files = []
    for sample in samples:
        if experiments.loc[sample, "condition"] == config["baseline_condition"]:
            baseline_samples.append(sample)
            baseline_files.append(experiments.loc[sample, "file"])
    return baseline_samples, baseline_files


def get_experiment_samples(experiments, samples):
    """Returns a list of experiment (i.e., non-baseline) sample and files"""
    experiment_samples = []
    experiment_files = []
    for sample in samples:
        if experiments.loc[sample, "condition"] != config["baseline_condition"]:
            experiment_samples.append(sample)
            experiment_files.append(experiments.loc[sample, "file"])
    return experiment_samples, experiment_files


def get_enrich2_input(wildcards):
    """Depending on whether or not filtering of zero counts is enabled,
    use the appropriate set of input counts for enrich2."""

    if remove_zeros:
        return expand(
            "results/{{experiment_name}}/processed_counts/removed_zeros/{samples}.tsv",
            samples=samples,
        )
    else:
        return expand(
            "results/{{experiment_name}}/processed_counts/{samples}.tsv",
            samples=samples,
        )


def get_input(wildcards):
    """Generate the input files for the dummy rule all.
    This is necessary to allow optional pipeline outputs."""

    input_list = []

    if experiment_samples:
        input_list.extend(
            expand(
                "results/{experiment_name}/rosace/{conditions}_scores.csv",
                experiment_name=config["experiment"],
                conditions=experimental_conditions,
            ) +
            expand(
                "results/{experiment_name}/rosace/rosace_installed.txt",
                experiment_name=config["experiment"]
            ),
        )
        if config["run_qc"]:
            input_list.extend(
                expand(
                    "stats/{experiment_name}/{experiment_name}_multiqc.html",
                    experiment_name=config["experiment"],
                )
            )
        if config["enrich2"]:
            input_list.extend(
                expand(
                    "results/{experiment_name}/enrich/tsv/{experiment_name}_exp/main_identifiers_scores.tsv",
                    experiment_name=config["experiment"],
                )
            )
    if config["baseline_condition"] and config["run_qc"]:
        input_list.extend(
            expand(
                "stats/{experiment_name}/{experiment_name}_baseline_multiqc.html",
                experiment_name=config["experiment"],
            )
        )

    return input_list


# Validate config and experiment files
validate(config, "../schemas/config.schema.yaml")

experiments = (
    pd.read_csv(config["experiment_file"], header=0)
    .dropna(how="all")
    .set_index("sample", drop=False, verify_integrity=True)
)

validate(experiments, "../schemas/experiments.schema.yaml")

# Set variables from config and experiment files
experiment = config["experiment"]
samples = experiments["sample"]
baseline_samples, baseline_files = get_baseline_samples(experiments, samples)
experiment_samples, experiment_files = get_experiment_samples(experiments, samples)
files = experiments["file"]
conditions = set(experiments["condition"])
experimental_conditions = conditions - set([config["baseline_condition"]])
reference_name = get_ref(config["reference"])
adapters_ref = pass_names(config["adapters"])
contaminants_ref = pass_names(config["contaminants"])
samtools_local = config["samtools_local"]

tiles = experiments["tile"].unique()
if len(tiles) > 1:
    config["tiled"] = True
else:
    config["tiled"] = False

# Do not remove unobserved variants if enrich will not be run
if config["enrich2"]:
    remove_zeros = config["remove_zeros"]
else:
    remove_zeros = False
