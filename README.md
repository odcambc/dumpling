# GATK-based Snakemake pipeline for deep mutational scanning experiments

This repository contains the Snakemake-based workflow for implementing
deep mutational scanning experiments used in the [Fraser](https://fraserlab.com/)
and [Coyote-Maestas](https://www.wcoyotelab.com/) labs.

Briefly, this conducts initial QC and mapping using BBTools, followed by the
AnalyzeSaturationMutagenesis GATK module to call variants in each replicate. After variant calling, the list of observed variants in each read is filtered
based on the list of designed variants, and the resulting counts are used to
infer the fitness of each variant using [Rosace](https://github.com/pimentellab/rosace)
and (optionally) [Enrich2](https://github.com/FowlerLab/Enrich2).

The pipeline is designed to be flexible and modular and should be amenable to use
with a variety of experimental designs. Please note several current [limitations](#limitations), however.

- [GATK-based Snakemake pipeline for deep mutational scanning experiments](#gatk-based-snakemake-pipeline-for-deep-mutational-scanning-experiments)
  - [Quick start](#quick-start)
  - [Installation](#installation)
    - [Install via GitHub](#install-via-github)
    - [Installing Rosace](#installing-rosace)
    - [Dependencies](#dependencies)
      - [Via conda (recommended)](#via-conda-recommended)
      - [Manually](#manually)
  - [Configuration](#configuration)
    - [Configuration files](#configuration-files)
    - [Working directory structure](#working-directory-structure)
  - [Usage](#usage)
    - [Running the pipeline](#running-the-pipeline)
    - [Output files](#output-files)
    - [Analyzing results](#analyzing-results)
      - [QC metrics](#qc-metrics)
      - [Data analysis](#data-analysis)
  - [Limitations](#limitations)
  - [Citations](#citations)
  - [License](#license)
  - [Contributing](#contributing)
  - [Getting help](#getting-help)

## Quick start

```bash
git clone https://github.com/odcambc/dumpling
cd dumpling
conda env create --file dumpling_env.yaml
conda activate dumpling_env
```

Note that, on ARM-based Macs, the conda environment may fail to install due to required packages not being available for that platform.
Assuming that [Rosetta](https://support.apple.com/en-us/102527) is installed, the environment can be installed using emulation with the following command:

```bash
CONDA_SUBDIR=osx-64 conda env create --file dumpling_env.yaml
```

If the environment installed and activated properly,
edit the configuration files in the `config` directory as needed. Then run the pipeline with:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```
## Installation

### Install via GitHub

Download or fork this repository and edit the configuration files as needed.

### Installing Rosace
This pipeline uses the [Rosace](https://github.com/pimentellab/rosace) scoring tool.
Rosace uses [CmdStanR](https://mc-stan.org/cmdstanr/) and R to infer scores.

Dumpling uses [renv](https://rstudio.github.io/renv/index.html) to handle R dependencies.
This pipeline also includes a minimal faculty to install Rosace automatically, but issues are
possible. This can be invoked by calling the `install_rosace` rule:

```snakemake --cores 8 install_rosace```

This tries to install renv, restore the renv environment, and install Rosace and CmdStanR. If this fails,
please try installing Rosace manually.

We recommend installing Rosace manually before running the pipeline, or at least 
verifying that the install script works. More details about manually installing Rosace are
available in the vignettes of the package and at the repository linked above.

### Dependencies

#### Via conda (recommended)
The simplest way to handle dependencies is with [Conda](https://conda.io/docs/) and the provided environment file.

```bash
conda env create --file dumpling_env.yaml
```

This will create a new environment named `dumpling` with all the dependencies installed. Then simply activate the environment and you're ready to go.

```bash
conda activate dumpling_env
```

#### Manually

The following are the dependencies required to run the pipeline:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [GATK](https://software.broadinstitute.org/gatk/)
* [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
* [Samtools](http://www.htslib.org/)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info/)
* [Enrich2](https://enrich2.readthedocs.io/en/latest/)
* [Rosace](https://github.com/pimentellab/rosace)

## Configuration

### Configuration files
The details of an experiment need to be specified in a configuration file that defines
parameters and an associated experiment file that details the experimental setup.

The configuration file is a YAML file: full details are included in the example file
`config/test_config.yaml` and in the schema file `schemas/config.schema.yaml`.

The experiment file is a CSV file that relates experimental conditions,
replicates, and time points to sequencing files: full details are included
in the config file and in the schema file `schemas/experiments.schema.yaml`.

Additionally, a reference fasta file is required for mapping. This should be
placed in the `references` directory, and the path to the file should be specified in the config file.

This pipeline also employs a processing step to standardize variant nomenclature
and remove any variants that are not designed or are likely errors. This
requires a CSV file containing the set of designed variants, including their
specific codon changes. This should be placed in the `config/designed_variants` directory,
and the path to the file should be specified in the config file. An example file
is included in `config/designed_variants/test_variants.csv`. This pipeline can generate
the variants CSV from the output set of oligos produced by the [DIMPLE](https://github.com/coywil26/DIMPLE)
library generation protocol: this can be enabled by including the path to the oligo CSV file in the config
file and setting `regenerate_variants` to `True` in the config.

### Working directory structure

The pipeline has the following directory structure:
```
├── workflow
│   ├── rules
│   ├── envs
│   ├── scripts
│   └── Snakefile
├── config
│   ├── test_config.yaml
│   ├── test_config.csv
│   ├── designed_variants
│   │   └── test_variants.csv
│   └── oligos
│       └── test_oligos.csv
├── logs
│   └── ...
├── references
│   └── test_ref.fasta
├── results
│   └── ...
├── schemas
│   ├── config.schema.yaml
│   └── experiments.schema.yaml
├── stats
│   └── ...
├── resources
│   ├── adapters.fa
│   ├── sequencing_artifacts.fa.gz
│   └── ...

```
## Usage
We normally use one instance of the pipeline for each experiment.
This allows for simpler tracking and reproducibility of individual experiments: for
a new dataset, fork the repo, edit the configuration files, and run the pipeline. This way,
a record of the exact configuration and environment can be saved. It is possible to run multiple
experiments in the same folder, but this is more difficult to reproduce.

### Running the pipeline

Once the dependencies have been installed (whether via conda or otherwise) the pipeline can be run with the following command:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 8
```

The maximum number of cores can be specified with the `--cores` flag. The `--software-deployment-method conda` flag 
tells Snakemake to use conda to create the environment specified within each rule.

### Output files

The pipeline generates a variety of output files. These are organized into the following directories:

* `benchmarks`: details of the runtime and process usage for each rule 
* `logs`: log files from each rule
* `results`: outputs from each rule (Note: many of these are intermediate files and are deleted by default).
* `stats`: various processing statistics from each rule
* `ref`: mapping target files generated by BBTools

These are ignored by git by default.

### Analyzing results
#### QC metrics
A variety of stats from tool outputs are provided in the `stats` directory. These are
aggregated using MultiQC. The aggregated reports contain:
* FastQC reports for raw reads (read counts, base quality, adapter content, etc.)
* BBTools reports
  * BBDuk reports for adapter trimming and contamination removal
  * BBMerge reports for merging paired-end reads
  * BBMap reports for mapping reads to the reference
* GATK AnalyzeSaturationMutagenesis reports for variant calling
* Reports for variant filtering

If a baseline condition is defined, a separate baseline report is also generated.

The files are saved as `stats/{experiment_name}_multiqc_report.html` and
`stats/{experiment_name}_baseline_multiqc_report.html` by default.

#### Data analysis

A starting analysis and plotting workflow is available in an associated 
repository: https://github.com/odcambc/dms_analysis_stub

## Limitations

We aim to regularly update this pipeline and continually expand
its functionality. However, there are currently several known limitations.

* The pipeline is currently designed for short-read sequencing. It does not support long-read PacBio or Nanopore sequencing.
* The pipeline is currently designed for direct sequencing. It does not support barcoded sequencing.
* The pipeline is currently designed for single-site variants (including varying-length indels, however). It largely does not support combinatorial variants.
* The designed variant generation step is currently optimized for DIMPLE libraries. Other protocols may require the user to generate the designed variants CSV themself.
* Rosace is designed for growth-based experiments. It is not optimized for FACS-seq experiments.

## Citations
This workflow is described in the following publication:

* [Rao et al., 2023](https://www.biorxiv.org/content/10.1101/2023.10.24.562292v1)

## License

This is licensed under the MIT license. See the LICENSE file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.wcoyotelab.com/members/).
