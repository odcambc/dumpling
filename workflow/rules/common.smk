from pathlib import Path
import logging


def get_file_from_sample(wildcards):
    """Maps from the sequencing output file names to the sample names defined in the experiment CSV.
    This is used in rule bbduk_trim_adapters"""

    data_dir = Path(config["data_dir"]).resolve()

    filename_match = experiments.loc[
        experiments["sample"] == wildcards.sample_prefix, "file"
    ]

    if filename_match.empty:
        raise ValueError(
            f"No matching file found for sample prefix: {wildcards.sample_prefix}"
        )

    if len(filename_match) > 1:
        raise ValueError(
            f"Multiple matching files found for sample prefix: {wildcards.sample_prefix}"
        )

    filename = filename_match.squeeze()

    # Find the actual file names for the given sample.
    R1, R2 = None, None
    for file in data_dir.glob(f"{filename}*"):
        if "_R1" in file.stem:
            R1 = file
        elif "_R2" in file.stem:
            R2 = file

    # Make sure both exist.
    if not R1 or not R2:
        raise FileNotFoundError(
            f"Could not find matching R1 and R2 files for prefix: {filename}"
        )

    return {"R1": str(R1), "R2": str(R2)}


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
            "results/{{experiment_name}}/processed_counts/enrich_format/{samples}.tsv",
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
            )
            + expand(
                "results/{experiment_name}/rosace/rosace_installed.txt",
                experiment_name=config["experiment"],
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


def validate_config(config):
    # Check ORF definition
    # etc etc...

    # Check reference file type
    reference_file_suffix = Path(config["reference"]).suffix
    if reference_file_suffix not in [
        ".fasta",
        ".fas",
        ".fa",
        ".fna",
        ".ffn",
        ".faa",
        ".mpfa",
        ".frn",
    ]:
        raise ValueError(
            f"Reference file {config['reference']} does not appear to be a valid FASTA file: {reference_file_suffix}"
        )

    # Check whether a correct set of input files are present
    if not Path(config["data_dir"]).exists():
        raise FileNotFoundError(f"Data directory {config['data_dir']} does not exist")

    if not reference_file.exists():
        raise FileNotFoundError(f"Reference file {config['reference']} does not exist")

    # Check if adapters is a list or a string, then check if each file exists
    adapter_paths = (
        config["adapters"]
        if isinstance(config["adapters"], list)
        else [config["adapters"]]
    )
    for adapter in adapter_paths:
        if not Path(adapter).exists():
            raise FileNotFoundError(f"Adapters file {adapter} does not exist")

    # Check if contaminants is a list or a string, then check if each file exists
    contaminants_paths = (
        config["contaminants"]
        if isinstance(config["contaminants"], list)
        else [config["contaminants"]]
    )
    for contaminant in contaminants_paths:
        if not Path(contaminant).exists():
            raise FileNotFoundError(f"Contaminants file {contaminant} does not exist")

    # Check for mode with variant filtering
    if not noprocess:
        # If we are regenerating the variants, then the variants file should not exist
        if config["regenerate_variants"]:
            if Path(variants_file).exists():
                print(
                    f"Variants file {variants_file} already exists, but regenerate_variants is set to true"
                )
            if not Path(oligo_file).exists():
                raise FileNotFoundError(
                    f"Oligo file {oligo_file} does not exist, but is required to make list of designed variants"
                )
        else:
            if not Path(variants_file).exists():
                raise FileNotFoundError(
                    f"Variants file {variants_file} does not exist, but is required to process counts"
                )

    return


# Validate config and experiment files
config.setdefault("bbtools_use_bgzip", True)
validate(config, "../schemas/config.schema.yaml")

experiments = (
    pd.read_csv(config["experiment_file"], header=0)
    .dropna(how="all")
    .set_index("sample", drop=False, verify_integrity=True)
)

validate(experiments, "../schemas/experiments.schema.yaml")

# Determine experiments, samples, and files by parsing config and experiment input
experiment = config["experiment"]
samples = experiments["sample"]
baseline_samples, baseline_files = get_baseline_samples(experiments, samples)
experiment_samples, experiment_files = get_experiment_samples(experiments, samples)
files = experiments["file"]
conditions = set(experiments["condition"])
experimental_conditions = conditions - set([config["baseline_condition"]])


# Load additional files
variants_file = config["variants_file"]
oligo_file = config["oligo_file"]

reference_file = Path(config["ref_dir"]) / config["reference"]
reference_name = Path(config["reference"]).stem

adapters_ref = pass_names(config["adapters"])
contaminants_ref = pass_names(config["contaminants"])

# Set configuration variables
samtools_local = config["samtools_local"]
noprocess = config["noprocess"]
bbtools_compression_flags = "" if config["bbtools_use_bgzip"] else "bgzip=f unbgzip=f "

# Set up tiled experiments
if "tile" not in experiments.columns:
    experiments["tile"] = 1

tiles = experiments["tile"].unique()

# If tiles are defined and they are all identical, then the experiment is not tiled
if len(tiles) > 1:
    config["tiled"] = True
else:
    config["tiled"] = False

# Do not remove unobserved variants if enrich will not be run
if config["enrich2"]:
    remove_zeros = config["remove_zeros"]
else:
    logging.warning("Enrich2 will not be run, so zero counts will not be removed.")
    remove_zeros = False

# Validate the configuration
validate_config(config)
