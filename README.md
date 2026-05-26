# GATK-based Snakemake pipeline for deep mutational scanning experiments

This repository contains the Snakemake-based workflow for implementing
deep mutational scanning experiments used in the [Fraser](https://fraserlab.com/)
and [Coyote-Maestas](https://www.wcoyotelab.com/) labs.

Briefly, this conducts initial QC and read processing with BBTools (adapter
trimming, contaminant filtering, and BBMerge error correction) and maps the
processed reads to the reference using either BBMap (default) or minimap2
(opt-in, configurable). The aligned reads are passed through GATK's
AnalyzeSaturationMutagenesis module to call variants in each replicate. After
variant calling, the list of observed variants in each read is filtered
based on the list of designed variants, and the resulting counts are used to
infer the fitness of each variant using [Rosace](https://github.com/pimentellab/rosace), [Lilace](https://github.com/pimentellab/lilace),
and (optionally) [Enrich2](https://github.com/FowlerLab/Enrich2).

The pipeline is designed to be flexible and modular and should be amenable to use
with a variety of experimental designs. Please note several current [limitations](#limitations), however.

- [GATK-based Snakemake pipeline for deep mutational scanning experiments](#gatk-based-snakemake-pipeline-for-deep-mutational-scanning-experiments)
  - [Quick start](#quick-start)
    - [Testing the pipeline with example data](#testing-the-pipeline-with-example-data)
  - [Installation](#installation)
    - [Install via GitHub](#install-via-github)
    - [Installing Rosace and Lilace](#installing-rosace-and-lilace)
      - [Issues installing Rosace and Lilace on OSX](#issues-installing-rosace-and-lilace-on-osx)
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

Installation for ARM-based Macs:

```bash
CONDA_SUBDIR=osx-64 conda env create --file dumpling_env.yaml
CONDA_SUBDIR=osx-64 conda env create --name enrich2 --file workflow/envs/enrich2.yaml

conda env create --platform osx-64 --name enrich2_arm64
```

You will also need to set the "samtools_local" variable in the config yaml to "true" to tell the pipeline to use this local version.

If the environment installed and activated properly,
edit the configuration files in the `config` directory as needed. Then run the pipeline with:

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 16
```

### Testing the pipeline with example data

To test the pipeline with example data and examine the output, you can
use the file provided in the [dumpling-example](https://github.com/odcambc/dumpling-example) repository. This repository contains a small dataset and configuration files that can be used to test the pipeline. To use it,
clone the repository and move the data directory into the dumpling directory, then create the environment and run the pipeline as above. The repository also
includes output files from running the pipeline on the example data that can
be used to compare results.

## Installation

### Install via GitHub

Download or fork this repository and edit the configuration files as needed.

### Installing Rosace, Lilace, and rosace-aa

This pipeline supports three scoring backends, all from the pimentellab group:
[Rosace](https://github.com/pimentellab/rosace), [Lilace](https://github.com/pimentellab/lilace), and
[rosace-aa](https://github.com/pimentellab/rosace-aa). All three use
[CmdStanR](https://mc-stan.org/cmdstanr/) and R to infer scores. rosace-aa is an extension of rosace that
decomposes the score into position + amino-acid substitution effects rather than a single per-variant scalar;
its score CSV layout matches rosace's so downstream tooling (e.g. `format_mavedb.py`) works unchanged.

Pick a backend via `scoring_backend: rosace | lilace | rosace_aa` in your config; `rosace` is the default.
Note that both `lilace` and `rosace_aa` require parsed variant metadata (wildtype/mutation/synonymous-control
columns) and are incompatible with `noprocess: true`.

Dumpling uses [renv](https://rstudio.github.io/renv/index.html) to handle R dependencies.
This pipeline also includes a minimal faculty to install each backend automatically, but issues are
possible. Invoke the relevant install rule:

```bash
snakemake --cores 8 install_rosace
snakemake --cores 8 install_lilace
snakemake --cores 8 install_rosace_aa
```

These try to install renv, restore the renv environment, and install the chosen backend with CmdStanR.
For `install_rosace_aa`, an additional `renv::install("pimentellab/rosace-aa@<sha>")` step pulls rosace-aa
from GitHub at a pinned SHA (the upstream repo has no tagged releases yet). If any install fails, please
try installing the package manually.

We recommend trying to install your chosen backend manually before running the pipeline, or at least
verifying that the install script works. More details about manual install are available in each
package's vignettes at the repository linked above.

#### Issues installing Rosace, Lilace, or rosace-aa on OSX

All three backends require a C++ and fortran compiler to install required dependencies.
R, by default, requires these to be installed in `/opt/gfortran`. User installs (via Homebrew, for example)
may not work. If you encounter an error compiling packages for the scoring backends, you may need to install
the gfortran compiler from R.

See https://cran.r-project.org/bin/macosx/tools/ for more details.

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

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [GATK](https://software.broadinstitute.org/gatk/)
- [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
- [minimap2](https://github.com/lh3/minimap2) (optional, opt-in alternative aligner — see [Aligner choice](#aligner-choice))
- [Samtools](http://www.htslib.org/)
- [pysam](https://github.com/pysam-developers/pysam)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](http://multiqc.info/)
- [Enrich2](https://enrich2.readthedocs.io/en/latest/)
- [Rosace](https://github.com/pimentellab/rosace)
- [Lilace](https://github.com/pimentellab/lilace)

BBTools compressed IO defaults to `pigz` (parallelized across each rule's threads —
typically saves 30-40 s/sample on 8 GB+ inputs vs single-threaded bgzip). Override via
`bbtools_compression: bgzip | pigz | none` in the config. `none` falls back to gzip and
is the right knob if your environment hangs in `bbduk.sh`, `bbmerge.sh`, or `bbmap.sh`.
The legacy `bbtools_use_bgzip: true|false` knob still works (with a deprecation warning)
and translates to `bgzip`/`none`.

## Configuration

### Configuration files

The details of an experiment need to be specified in a configuration file that defines
parameters and an associated experiment file that details the experimental setup.

The configuration file is a YAML file: full details are included in the example file
`config/example.yaml` and in the schema file `workflow/schemas/config.schema.yaml`.

The experiment file is a CSV file that relates experimental conditions,
replicates, and time points to sequencing files: full details are included
in the config file and in the schema file `workflow/schemas/experiments.schema.yaml`.

Additionally, a reference fasta file is required for mapping. This should be
placed in the `references` directory, and the path to the file should be specified in the config file.

This pipeline also employs a processing step to standardize variant nomenclature
and remove any variants that are not designed or are likely errors. This
requires a CSV file containing the set of designed variants, including their
specific codon changes. This should be placed in the `config/designed_variants` directory,
and the path to the file should be specified in the config file. An example file
is included in `config/designed_variants/example_variants.csv`. This pipeline can generate
the variants CSV from the output set of oligos produced by the [DIMPLE](https://github.com/coywil26/DIMPLE)
library generation protocol: this can be enabled by including the path to the oligo CSV file in the config
file and setting `regenerate_variants` to `True` in the config.

### Aligner choice

The pipeline supports two alternative aligners for mapping reads to the
reference, selected via the `aligner` key in the config:

```yaml
aligner: bbmap      # default — current behavior
# aligner: minimap2 # opt-in
```

| Aligner | Wall (example fixture) | Peak RSS | Mapping-stage QC artifacts |
|---|---|---|---|
| `bbmap` (default) | ~37 s/sample | ~8.5 GB | BBMap-format histograms (`_map.covstats`, `_map.basecov`, `_map.ehist`, `_map.indelhist`, `_map.mhist`, `_map.idhist`, `_map.bincov`) |
| `minimap2` | ~5 s/sample | ~450 MB | samtools-format outputs (`_samtools_stats`, `_samtools_flagstat`) |

Both produce biologically equivalent variant counts (Rosace score Pearson
r > 0.997 between the two on the example fixture). At production scale
(multi-GB FASTQs) the wall advantage of minimap2 shrinks to ~2-3× as BBMap's
fixed index-loading overhead amortizes; the RSS advantage holds.

Per-position coverage stats are not currently produced under `minimap2`
(`samtools coverage` is too expensive on deeply-covered DMS data); a
mosdepth-based replacement is planned. All other MultiQC sections work
identically under both aligners.

### Working directory structure

The pipeline has the following directory structure:

```
├── workflow
│   ├── rules
│   ├── envs
│   ├── scripts
│   ├── schemas
│   │   ├── config.schema.yaml
│   │   └── experiments.schema.yaml
│   └── Snakefile
├── config
│   ├── example.yaml
│   ├── example.csv
│   ├── multiqc_config.yaml
│   ├── designed_variants
│   │   └── example_variants.csv
│   └── oligos
│       └── (optional DIMPLE oligo CSVs)
├── logs
│   └── ...
├── references
│   └── example_ref.fasta
├── results
│   └── ...
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

- `benchmarks`: details of the runtime and process usage for each rule
- `logs`: log files from each rule
- `results`: outputs from each rule (Note: many of these are intermediate files and are deleted by default).
- `stats`: various processing statistics from each rule
- `ref`: normalized reference, sequence dictionary, and aligner indexes. Reference-derived
  artifacts here persist across pipeline runs — keyed on a content hash of the reference, so
  unchanged references skip the index rebuild on repeat runs.

These are ignored by git by default.

Sample-derived intermediates (trimmed/cleaned FASTQs under `results/{experiment}/` and the
mapped BAMs) are marked `temp()` and deleted as soon as their downstream consumers finish.
To retain them for debugging — e.g. inspecting a `bbduk`-trimmed FASTQ or running
`samtools view` on a mapped BAM — pass `--notemp` (alias: `--no-temp`, `--nt`) to the
snakemake invocation:

```bash
snakemake -s workflow/Snakefile --configfile config/my_experiment.yaml --cores 16 --notemp
```

`--notemp` is a Snakemake built-in; no config knob needed.

### Analyzing results

#### QC metrics

A variety of stats from tool outputs are provided in the `stats` directory. These are
aggregated using MultiQC. The aggregated reports contain:
- FastQC reports for raw reads (read counts, base quality, adapter content, etc.)
- BBTools reports
  - BBDuk reports for adapter trimming and contamination removal
  - BBMerge reports for merging paired-end reads
  - BBMap reports for mapping reads to the reference *(when `aligner: bbmap`)*
- samtools stats and flagstat reports *(when `aligner: minimap2`)*
- GATK AnalyzeSaturationMutagenesis reports for variant calling
- Reports for variant filtering

If a baseline condition is defined, a separate baseline report is also generated.

The files are saved as `stats/{experiment_name}_multiqc_report.html` and
`stats/{experiment_name}_baseline_multiqc_report.html` by default.

#### Data analysis

A starting analysis and plotting workflow is available in an associated
repository: <https://github.com/odcambc/dms_analysis_stub>

## Limitations

We aim to regularly update this pipeline and continually expand
its functionality. However, there are currently several known limitations.

- The pipeline is currently designed for short-read sequencing. It does not support long-read PacBio or Nanopore sequencing.
- The pipeline is currently designed for direct sequencing. It does not support barcoded sequencing.
- The pipeline is currently designed for single-site variants (including varying-length indels, however). It largely does not support combinatorial variants.
- The designed variant generation step is currently optimized for DIMPLE libraries. Other protocols may require the user to generate the designed variants CSV themself.
- This pipeline may not work properly if the data is in a cloud server (i.e., a Box drive) or other non-standard file system.
- This pipeline currently only accepts fastq.gz files. It does not accept fastq files.

## Citations

This workflow, along with Rosace, is described in the following publication:

- Preprint: [Rao et al., 2023](https://www.biorxiv.org/content/10.1101/2023.10.24.562292v1)
- Published: [Rao et al., 2024](https://doi.org/10.1186/s13059-024-03279-7)

The Lilace FACS-based model is described in:

- [Freudenberg et al., 2026](https://doi.org/10.1186/s13059-026-03934-1)

## License

This is licensed under the MIT license. See the LICENSE file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.wcoyotelab.com/members/).
