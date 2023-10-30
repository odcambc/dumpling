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

# Validate config and experiment files
validate(config, "../schemas/config.schema.yaml")

experiments = pd.read_csv(config["experiment_file"], header=0).set_index(
    "sample", drop=False, verify_integrity=True
)

validate(experiments, "../schemas/experiments.schema.yaml")

# Set variables from config and experiment files
experiment = config["experiment"]
samples = experiments["sample"]
files = experiments["file"]
conditions = set(experiments["condition"])
reference_name = get_ref(config["reference"])
adapters_ref = pass_names(config["adapters"])
contaminants_ref = pass_names(config["contaminants"])