# A Snakemake pipeline for deep mutational scanning experiments

This repository contains Dumpling, the Snakemake-based workflow for implementing
deep mutational scanning experiments used in the [Fraser](https://fraserlab.com/)
and [Coyote-Maestas](https://www.wcoyotelab.com/) labs.

## What is dumpling

Dumpling is an end-to-end pipeline for analyzing and scoring deep mutational scanning experiments.

Dumpling will perform sequence QC, variant calling, filtering, and scoring. We made Dumpling to complement [DIMPLE](https://github.com/coywil26/DIMPLE), our library generation platform, but it should be more generally useful for other types of DMS libraries. DIMPLE allows insertion and deletion variants in addition to substitutions, which dumpling was explicitly created to support.

Because DIMPLE libraries generate libraries of specifically-designed variants, dumpling by default filters identified variants based on an expected set of inputs: this can be easily disabled to support random mutagenesis.

Dumpling supports the following sequencing strategies:

* Tagmentation
* Tiled amplicons

Dumpling supports a number of experimental designs:

* Pooled growth and competitions
* FACS screens
* Binding and competition

We are interested and happy to discuss accommodating new capabilities - drop us a line if you're interested!

Dumpling does _not_ currently support barcoded or UMI libraries and has limited support for
combinatorial or multi-mutation libraries. It also does not currently support insertion scanning libraries or fusion libraries. These are planned for future releases.

Dumpling was also designed to complement the scoring tools we developed in the lab, and so provides end-to-end scoring for
[Rosace](https://github.com/pimentellab/rosace),
[Rosace-AA](https://github.com/pimentellab/rosace-aa), and
[Lilace](https://github.com/pimentellab/lilace). [Enrich2](https://github.com/FowlerLab/Enrich2) scores are also included.

Dumpling also supports [Cosmos](https://github.com/pimentellab/cosmos) for multi-phenotype causal modelling after individual experiment scoring.

A [configuration generation tool](https://dumpling.odcambc.com) is available too!

## Quick start

A prebuilt container image is published to [GitHub Container Registry](https://github.com/odcambc/dumpling/pkgs/container/dumpling) which bundles Snakemake and the full dumpling toolchain — BBTools, GATK, minimap2, samtools, FastQC, MultiQC, R 4.5 with the Rosace/Lilace/Rosace-AA scorers, **CmdStan pre-compiled**, and Enrich2. This is the fastest path to a working pipeline because it skips conda environment setup and tool installation, which are most likely to hit user environment issues.

The same image works for both Docker (local dev, CI) and Apptainer/Singularity (HPC), since Apptainer transparently pulls and converts OCI images.

### 1. Clone and write your configuration

​```bash
git clone https://github.com/odcambc/dumpling
cd dumpling
​```

dumpling needs three inputs, all referenced from your config YAML (see
[config/README.md](config/README.md) for the full reference):

* a **config file** — copy `config/example.yaml` to `config/my_experiment.yaml` and edit, or use the [configuration generator](https://dumpling.odcambc.com)
* an **experiment CSV** — maps samples/replicates/timepoints to FASTQ files
* a **reference FASTA** under `references/`

Set `data_dir: data` in your config — that's where your reads get mounted below.

### 2. Run it with Docker

Your FASTQs usually live outside the repo. Mount the repo at `/workdir` (the
container's working directory) and your reads at `/workdir/data` (matching
`data_dir`).

​```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd):/workdir" \
  -v "/path/to/your/fastqs:/workdir/data" \
  ghcr.io/odcambc/dumpling:latest \
  snakemake --configfile config/my_experiment.yaml --cores 16
​```

Outputs (`results/`, `stats/`, `logs/`) land back in the repo via the `$(pwd)`
mount; `--user` keeps them owned by you rather than root.

## Installation

### With Conda (locally)

If you'd rather manage dependencies yourself, the conda-based install is below.

```bash
git clone https://github.com/odcambc/dumpling
cd dumpling
conda env create --file dumpling_env.yaml
conda activate dumpling_env
```

See [Troubleshooting](#troubleshooting) for issues running dumpling on ARM Macs or installing Rosace.

### Dependencies

The following are the dependencies required to run the pipeline:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [GATK](https://software.broadinstitute.org/gatk/)
* [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
* [minimap2](https://github.com/lh3/minimap2) (optional, opt-in alternative aligner — see [Aligner choice](#aligner-choice))
* [Samtools](http://www.htslib.org/)
* [pysam](https://github.com/pysam-developers/pysam)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info/)
* [Enrich2](https://enrich2.readthedocs.io/en/latest/)
* [Rosace](https://github.com/pimentellab/rosace)
* [Rosace-AA](https://github.com/pimentellab/rosace-aa)
* [Lilace](https://github.com/pimentellab/lilace)

## Configuration

Try the online [configuration generation tool](https://dumpling.odcambc.com)!

### Input files

dumpling requires three things: a configuration file, an experiment/sample definition file, and a reference sequence fasta.

To use variant filtering, it also needs a list of designed variants. This can be generated automatically from
a DIMPLE-generated list of oligos, which can also be provided.

Full details of configuration files are in the [configuration README](config/README.md).
The `config/example.yaml` file is populated with reasonable defaults if you just want to get going.

### Aligner choice

The pipeline supports two alternative aligners for mapping reads to the
reference, selected via the `aligner` key in the config:

```yaml
aligner: bbmap      # default — current behavior
# aligner: minimap2 # opt-in
```

minimap2 is significantly faster (up to an order of magnitude) in our experience but bbmap returns slightly more detailed QC metrics and is preserved as default. In testing the results are > 0.997 identical, which likely is simply non-determinism of multi-threaded mapping.

## Usage

We normally use one instance of the pipeline for each experiment.
This allows for simpler tracking and reproducibility of individual experiments: for
a new dataset, fork the repo, edit the configuration files, and run the pipeline. This way,
a record of the exact configuration and environment can be saved. It is possible to run multiple
experiments in the same folder, but this is more difficult to reproduce.

### Running the pipeline

Once the dependencies have been installed (whether via conda or otherwise) the pipeline can be run with the following command:

```bash
snakemake --configfile config/my_experiment.yaml --software-deployment-method conda --cores 8
```

The maximum number of cores can be specified with the `--cores` flag. The `--software-deployment-method conda` flag
tells Snakemake to use conda to create the environment specified within each rule.

For a local run that also respects per-rule memory budgets, use the bundled local profile:

```bash
snakemake --profile workflow/profiles/default --configfile config/my_experiment.yaml --cores 16
```

### With Apptainer/Singularity (HPC)

​```bash
snakemake --use-singularity --configfile config/my_experiment.yaml --cores 16
​```

Apptainer auto-binds your working directory, so keep `data_dir` reachable from
there (under the repo, or add `--singularity-args "--bind /path/to/fastqs:/workdir/data"`).
The `containerized:` directive in `workflow/Snakefile` points at the GHCR image,
so `--use-singularity` picks it up automatically. You do **not** need `--use-conda` —
the container bundles every tool in a single environment.

Available tags: `:latest` (most recent release), `:vX.Y.Z` (pinned release).

### Running on a cluster

Note: the following is a work in progress and may not work for your cluster or environment.

dumpling ships a SLURM profile (`workflow/profiles/slurm`). Snakemake runs on the
login node and submits one `sbatch` job per rule instance; each job executes
inside the prebuilt container on the compute nodes. Every heavy rule declares `threads` and `resources: mem_mb` (tunable via the `mem_*` config knobs), which Snakemake translates into `--cpus-per-task` / `--mem` / `--time`.

1. On the login node, create the thin submitting environment (Snakemake + the
   SLURM executor plugin only):

   ```bash
   conda env create -f cluster_env.yaml
   conda activate dumpling_cluster
   ```

2. Set your site's SLURM account and partition by uncommenting `slurm_account`
   / `slurm_partition` under `default-resources` in `workflow/profiles/slurm/config.yaml`
   (these are the only values that can't be defaulted). Adjust `mem_*` budgets in
   your config if the defaults don't fit your data — read `benchmarks/{experiment}/`
   `max_rss` from a real run to size them.

3. Launch:

   ```bash
   snakemake --profile workflow/profiles/slurm --configfile config/my_experiment.yaml --software-deployment-method apptainer
   ```

To validate the profile without submitting anything (e.g. to check resources and
account settings resolve), add `--dry-run`.

### Output files

The pipeline generates a variety of output files. These are organized into the following directories:

* `benchmarks`: details of the runtime and process usage for each rule
* `logs`: log files from each rule
* `results`: outputs from each rule (Note: many of these are intermediate files and are deleted by default).
* `stats`: various processing statistics from each rule
* `ref`: normalized reference, sequence dictionary, and aligner indexes. Reference-derived
  artifacts here persist across pipeline runs — keyed on a content hash of the reference, so
  unchanged references skip the index rebuild on repeat runs.

These are ignored by git by default.

Sample-derived intermediates (trimmed/cleaned FASTQs under `results/{experiment}/` and the
mapped BAMs) are marked `temp()` and deleted as soon as their downstream consumers finish.
To retain them for debugging — e.g. inspecting a `bbduk`-trimmed FASTQ or running
`samtools view` on a mapped BAM — pass `--notemp` (alias: `--no-temp`, `--nt`) to the
snakemake invocation:

```bash
snakemake --configfile config/my_experiment.yaml --cores 16 --notemp
```

`--notemp` is a Snakemake built-in; no config knob needed.

### Analyzing results

#### QC metrics

A variety of stats from tool outputs are provided in the `stats` directory. These are
aggregated using MultiQC. The aggregated reports contain:

* FastQC reports for raw reads (read counts, base quality, adapter content, etc.)
* BBTools reports
  * BBDuk reports for adapter trimming and contamination removal
  * BBMerge reports for merging paired-end reads
  * BBMap reports for mapping reads to the reference _(when `aligner: bbmap`)_
* samtools stats and flagstat reports _(when `aligner: minimap2`)_
* GATK AnalyzeSaturationMutagenesis reports for variant calling
* Reports for variant filtering

If a baseline condition is defined, a separate baseline report is also generated.

The files are saved as `stats/{experiment_name}_multiqc_report.html` and
`stats/{experiment_name}_baseline_multiqc_report.html` by default.

## Troubleshooting

### Using Conda on ARM Macs

Note that, on ARM-based Macs, the conda environment may fail to install due to required packages not being available for that platform. Compatibility is a moving target however, and this may not be accurate.

Assuming that [Rosetta](https://support.apple.com/en-us/102527) is installed, the environment can be installed using emulation with the following command:

```bash
CONDA_SUBDIR=osx-64 conda env create --file dumpling_env.yaml
CONDA_SUBDIR=osx-64 conda env create --name enrich2 --file workflow/envs/enrich2.yaml

conda env create --platform osx-64 --name enrich2_arm64
```

You will also need to set the "samtools_local" variable in the config yaml to "true" to tell the pipeline to use this local version.

If the environment installed and activated properly,
edit the configuration files in the `config` directory as needed, then run the pipeline.

### Installing Rosace, Lilace, and Rosace-AA

This pipeline supports three scoring backends, all from the pimentellab group:
[Rosace](https://github.com/pimentellab/rosace), [Lilace](https://github.com/pimentellab/lilace), and
[Rosace-AA](https://github.com/pimentellab/rosace-aa). All three use
[CmdStanR](https://mc-stan.org/cmdstanr/) and R to infer scores. Rosace-AA is an extension of rosace that
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
For `install_rosace_aa`, an additional `renv::install("pimentellab/rosace-aa@<sha>")` step pulls Rosace-AA
from GitHub at a pinned SHA (the upstream repo has no tagged releases yet). If any install fails, please
try installing the package manually.

We recommend trying to install your chosen backend manually before running the pipeline, or at least
verifying that the install script works. More details about manual install are available in each
package's vignettes at the repository linked above.

#### Issues installing Rosace, Lilace, or Rosace-AA on OSX

All three backends require a C++ and fortran compiler to install required dependencies.
R, by default, requires these to be installed in `/opt/gfortran`. User installs (via Homebrew, for example)
may not work. If you encounter an error compiling packages for the scoring backends, you may need to install
the gfortran compiler from R.

See <https://cran.r-project.org/bin/macosx/tools/> for more details.

## Reference

### Working directory structure

The pipeline has the following directory structure:

```
├── workflow
│   ├── rules
│   │   └── scripts
│   ├── envs
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

### Mapping tuning

BBTools compressed IO defaults to `pigz` (parallelized across each rule's threads —
typically saves 30-40 s/sample on 8 GB+ inputs vs single-threaded bgzip). Override via
`bbtools_compression: bgzip | pigz | none` in the config. `none` falls back to gzip and
is the right knob if your environment hangs in `bbduk.sh`, `bbmerge.sh`, or `bbmap.sh`.
The legacy `bbtools_use_bgzip: true|false` knob still works (with a deprecation warning)
and translates to `bgzip`/`none`.

## Limitations

We aim to regularly update this pipeline and continually expand
its functionality. However, there are currently several known limitations.

* The pipeline is currently designed for short-read sequencing. It does not support long-read PacBio or Nanopore sequencing.
* The pipeline is currently designed for direct sequencing. It does not support barcoded sequencing.
* The pipeline is currently designed for single-site variants (including varying-length indels, however). It largely does not support combinatorial variants.
* The designed variant generation step is currently optimized for DIMPLE libraries. Other protocols may require the user to generate the designed variants CSV themself.
* This pipeline may not work properly if the data is in a cloud server (i.e., a Box drive) or other non-standard file system.
* This pipeline currently only accepts fastq.gz files. It does not accept fastq files.

## Citations

This workflow, along with Rosace, is described in the following publication:

* Preprint: [Rao et al., 2023](https://www.biorxiv.org/content/10.1101/2023.10.24.562292v1)
* Published: [Rao et al., 2024](https://doi.org/10.1186/s13059-024-03279-7)

The Rosace-AA extension to Rosace is described in:

* Preprint: [Rao et al., 2025](https://www.biorxiv.org/content/10.1101/2025.01.09.632281v1)
* Published: [Rao et al., 2025](https://doi.org/10.1093/bioadv/vbaf218)

The Lilace FACS-based model is described in:

* Preprint: [Freudenberg et al., 2025](https://www.biorxiv.org/content/10.1101/2025.06.24.661380v1)
* Published: [Freudenberg et al., 2026](https://doi.org/10.1186/s13059-026-03934-1)

The Cosmos causal model is described in:

* Preprint: [Rao et al., 2025](https://www.biorxiv.org/content/10.1101/2025.08.01.667517v2)

The Enrich 2 model and tool is described in:

* Preprint: [Rubin et al., 2016](https://www.biorxiv.org/content/10.1101/075150v1.abstract)
* Published: [Rubin et al., 2017](https://doi.org/10.1186/s13059-017-1272-5)

## Other comparable tools

* [Enrich2](https://github.com/fowlerlab/enrich2)
* [DiMSum](https://github.com/lehner-lab/DiMSum)
* [gyōza](https://github.com/durr1602/gyoza)
* [mutscan](https://github.com/fmicompbio/mutscan)
* [ACIDES](https://github.com/nemoto-lab/ACIDES)

See [Çubuk et al., 2025](https://doi.org/10.1038/s44320-025-00137-x) for a good review of these and other DMS scoring approaches.

## License

This is licensed under the MIT license. See the LICENSE file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.waymentsteelelab.org).
