from pathlib import Path
import logging
import re

# Regex patterns for identifying paired-end read files.
# Matches common conventions: _R1_001, _R1, _1 with .fastq.gz, .fq.gz, .fastq, .fq extensions
_FASTQ_R1_PATTERN = re.compile(r"[._](?:R1|1)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")
_FASTQ_R2_PATTERN = re.compile(r"[._](?:R2|2)(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$")


def resolve_fastq_pair(data_dir, filename):
    """Resolve a file prefix to a paired-end R1/R2 fastq pair.

    Handles common naming conventions:
      - Illumina standard: {prefix}_R1_001.fastq.gz
      - Simplified: {prefix}_R1.fastq.gz
      - Numeric: {prefix}_1.fastq.gz
      - Any of the above with .fq.gz, .fastq, or .fq extensions
    """
    data_dir = Path(data_dir).resolve()
    R1, R2 = None, None

    for file in data_dir.glob(f"{filename}*"):
        name = file.name
        if _FASTQ_R1_PATTERN.search(name):
            if R1 is not None:
                raise ValueError(
                    f"Multiple R1 files found for prefix '{filename}': {R1.name} and {name}"
                )
            R1 = file
        elif _FASTQ_R2_PATTERN.search(name):
            if R2 is not None:
                raise ValueError(
                    f"Multiple R2 files found for prefix '{filename}': {R2.name} and {name}"
                )
            R2 = file

    if not R1 or not R2:
        raise FileNotFoundError(
            f"Could not find matching R1 and R2 fastq files for prefix '{filename}' in {data_dir}. "
            f"Expected files matching patterns like {filename}_R1_001.fastq.gz, {filename}_R1.fq.gz, "
            f"{filename}_1.fastq, etc."
        )

    return R1, R2


def get_file_from_sample(wildcards):
    """Maps from the sequencing output file names to the sample names defined in the experiment CSV.
    This is used in rule bbduk_trim_adapters"""

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
    R1, R2 = resolve_fastq_pair(config["data_dir"], filename)

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

# Resolve fastq file paths for all samples at parse time
fastq_map = {}
for file_prefix in files:
    try:
        r1, r2 = resolve_fastq_pair(config["data_dir"], file_prefix)
        fastq_map[file_prefix] = {"R1": r1, "R2": r2}
    except (FileNotFoundError, ValueError) as e:
        logging.warning(f"Could not resolve fastq pair for '{file_prefix}': {e}")

# Build list of raw fastq paths for QC rules
raw_fastq_paths = []
for file_prefix in files:
    if file_prefix in fastq_map:
        raw_fastq_paths.append(str(fastq_map[file_prefix]["R1"]))
        raw_fastq_paths.append(str(fastq_map[file_prefix]["R2"]))


def get_fastqc_name(fastq_path):
    """Derive the FastQC output name from a fastq file path.
    FastQC strips .gz and .fastq/.fq to form the base name."""
    name = Path(fastq_path).name
    name = re.sub(r"\.gz$", "", name)
    name = re.sub(r"\.(?:fastq|fq)$", "", name)
    return name


# Map file prefixes to their FastQC output base names
fastqc_names = {}
for file_prefix, paths in fastq_map.items():
    fastqc_names[file_prefix] = {
        "R1": get_fastqc_name(paths["R1"]),
        "R2": get_fastqc_name(paths["R2"]),
    }

# Reverse map: FastQC base name -> actual fastq file path (for fastqc rule input)
fastqc_input_map = {}
for file_prefix, paths in fastq_map.items():
    fastqc_input_map[get_fastqc_name(paths["R1"])] = str(paths["R1"])
    fastqc_input_map[get_fastqc_name(paths["R2"])] = str(paths["R2"])

# Validate the configuration
validate_config(config)
