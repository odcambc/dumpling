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

def get_input(wildcards):
    """Generate the input files for the dummy rule all.
    This is necessary to allow optional pipeline outputs."""

    input_list = []

    if experiment_samples:
        input_list.extend(
            expand(
                "results/{experiment_name}/rosace/{conditions}_scores.tsv",
                experiment_name=config["experiment"],
                conditions=conditions,
            )
        )
        input_list.extend(
            expand(
                "stats/{experiment_name}/{experiment_name}_multiqc.html",
                experiment_name=config["experiment"],
            )
        )
    if config["baseline_condition"]:
        input_list.extend(
            expand(
                "stats/{experiment_name}/{experiment_name}_baseline_multiqc.html",
                experiment_name=config["experiment"],
            )
        )

    return input_list

# Validate config and experiment files
validate(config, "../schemas/config.schema.yaml")

experiments = pd.read_csv(config["experiment_file"], header=0).dropna(how = 'all').set_index(
    "sample", drop=False, verify_integrity=True
)

validate(experiments, "../schemas/experiments.schema.yaml")

# Set variables from config and experiment files
experiment = config["experiment"]
samples = experiments["sample"]
baseline_samples, baseline_files = get_baseline_samples(experiments, samples)
experiment_samples, experiment_files = get_experiment_samples(experiments, samples)
files = experiments["file"]
conditions = set(experiments["condition"])
reference_name = get_ref(config["reference"])
adapters_ref = pass_names(config["adapters"])
contaminants_ref = pass_names(config["contaminants"])