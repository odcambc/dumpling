# GATK based snakemake pipeline for deep mutational scanning experiments

This repository contains the snakemake-based workflow for implementing
deep mutational scanning experiments used in the Fraser and Coyote-Maestas labs.

Briefly, this conducts initial QC and mapping using bbtools, followed by the
AnalyzeSaturationMutagenesis GATK module to call variants in each replicate.

## Quick start

```bash
git clone https://github.com/odcambc/dumpling
conda create --name dumpling --file spec-file.txt
conda activate dumpling
```
Edit the configuration files in the `config` directory as needed. Then run the pipeline with:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 8
```

## Installation

### Install via github

Download or fork this repository and edit the configuration files as needed.

### Dependencies

#### Via conda (recommended)
The simplest way to handle dependencies is with [Conda](https://conda.io/docs/) and the provided environment file.

```bash
conda create --name dumpling --file spec-file.txt
```

This will create a new environment named `dumpling` with all the dependencies installed. Then simply activate the environment and you're ready to go.

```bash
conda activate dumpling
```

#### Manually

The following are the dependencies required to run the pipeline:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [GATK](https://software.broadinstitute.org/gatk/)
* [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/)
* [samtools](http://www.htslib.org/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [multiqc](http://multiqc.info/)
* [Enrich2](https://enrich2.readthedocs.io/en/latest/)

## Configuration

### Configuration files
The details of an experiment need to be specified in a configuration file which defines parameters and an associated experiment file which details
the experimental setup.

The configuration file is a YAML file: full details are included in the example file `config/test_config.yaml` and in the schema file `schemas/config.schema.yaml`.

The experiment file is a CSV file which relates experimental conditions, replicates, and timepoints to sequencing files: full details are included in the config file and in the schema file `schemas/experiments.schema.yaml`.

Additionally, a reference fasta file is required for mapping. This should be placed in the `references` directory, and the path to the file should be specified in the config file.

This pipeline also employs a processing step to standardize the variant
nomenclature as well as remove any variants that are not designed or are likely errors. This requires a csv file containing the set of designed variants, including their specific codon changes. This should be placed in the `config/designed_variants` directory, and the path to the file should be specified in the config file. An example file is included in `config/designed_variants/test_variants.csv`. This pipeline also includes a faculty
to generate the variants csv from a set of oligos used to generate the library: this can be enabled by including the path to the oligo csv file in the config file, and setting `regenerate_variants` to `True` in the config.

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

### Running the pipeline

Once the dependencies have been installed (whether via conda or otherwise) the pipeline can be run with the following command:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 8
```

The maximum number of cores can be specified with the `--cores` flag. The `--software-deployment-method conda` flag tells snakemake to use conda to create the environment specified within each rule.

### Evaluating statistics

The pipeline creates a set of QC metrics for each sample, going from
raw reads (with FastQC) through to variant calling and scoring. These
are aggregated into a single report using MultiQC, which is saved in
`stats/{experiment_name}_multiqc_report.html` upon completion.
