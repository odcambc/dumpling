# GATK based snakemake pipeline for deep mutational scanning experiments

This repository contains the snakemake-based workflow for implementing
deep mutational scanning experiments used in the Fraser and Coyote-Maestes labs.

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
snakemake -s workflow/Snakefile --use-conda --cores 8
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

### Configuration files

## Usage

### Running the pipeline

Once the dependencies have been installed (whether via conda or otherwise) the pipeline can be run with the following command:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 8
```

The maximum number of cores can be specified with the `--cores` flag. The `--use-conda` flag tells snakemake to use conda to create the environment specified within each rule.

### Evaluating statistics

The pipeline creates a set of QC metrics for each sample, going from
raw reads (with FastQC) through to variant calling and scoring. These
are aggregated into a single report using MultiQC, which is saved in
`stats/{experiment_name}_multiqc_report.html` upon completion.
