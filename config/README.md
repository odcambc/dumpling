# Dumpling configuration

A dumpling run is defined by a **config YAML** plus the files it points to. Copy
[`example.yaml`](example.yaml) to `config/<your_experiment>.yaml`, edit, and pass
it with `--configfile`. Every key is validated against
[`../workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml);
a [configuration generator](https://dumpling.odcambc.com) is also available.

The [`example.yaml`](example.yaml) should have sensible defaults: a full list of options is listed below.

## Inputs: paths and files

| Key               | Meaning                                                            |
| ----------------- | ------------------------------------------------------------------ |
| `experiment`      | Unique experiment name; names output dirs/files.                   |
| `data_dir`        | Directory holding **all** raw reads (flat, not in subfolders).     |
| `ref_dir`         | Directory holding reference FASTA.                                 |
| `reference`       | Reference FASTA filename (nucleotide).                             |
| `experiment_file` | Experiment CSV (see below).                                        |
| `variants_file`   | Designed-variants CSV (see below).                                 |
| `oligo_file`      | *(optional)* DIMPLE oligo CSV, used to regenerate `variants_file`. |
| `orf`             | ORF coordinate range within the reference, e.g. `"141-1568"`.      |

### Experiment CSV

This file defines the experimental organization of the data.

Columns: `sample` (unique), `condition`, `replicate`, `time` (or `bin` for FACS),
`tile` number (for tiled amplicon sequencing), `file` (filename prefix, i.e. without
`_R1_001.fastq.gz`/`_R2_001.fastq.gz`). For `run_cosmos`, exactly two conditions
additionally set a `phenotype` column to `1` and `2` (the two sequential phenotypes
cosmos models). Validated against
[`../workflow/schemas/experiments.schema.yaml`](../workflow/schemas/experiments.schema.yaml).

### Reference FASTA

A nucleotide FASTA covering the expected mapped region, placed in `ref_dir`
(default `references/`) and named via the `reference` key. Set the ORF coordinates
within this file via `orf`: these are **1-referenced** (identical to SnapGene
numbering) and are required for correct variant calling.

### Designed-variants CSV

dumpling standardizes variant nomenclature and drops variants that aren't
designed (or are likely errors). Columns: `count` (init 0), `pos`,
`mutation_type` (`S`/`M`/`D`/`I`/`X`), `name` (e.g. `A123T`), `codon`,
`mutation` (subtype incl. indel length, e.g. `D_3`), `length`, `hgvs`.

Generate it from a DIMPLE `oligo_file` by setting `regenerate_variants: true`.
To skip designed-variant filtering entirely (e.g. random mutagenesis), set
`noprocess: true`.

Note that even if filtering is enabled, the filtered variant counts will be
saved as a "rejected" counts file.

## What to run

| Key                   | Default  | Effect                                                                                                                                                                                  |
| --------------------- | -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scoring_backend`     | `rosace` | Scoring model: `rosace`, `lilace`, or `rosace_aa`. `lilace`/`rosace_aa` need parsed variant metadata → incompatible with `noprocess: true`.                                             |
| `lilace_seed`         | `null`   | Random seed for the lilace backend; `null` runs nondeterministically. No effect unless `scoring_backend: lilace`.                                                                       |
| `enrich2`             | `true`   | Also run Enrich2 alongside the backend.                                                                                                                                                 |
| `keep_enrich_h5`      | `false`  | Keep Enrich2 `.h5` stores (large; otherwise `temp()`). No effect unless `enrich2`.                                                                                                      |
| `deposit_to_mavedb`   | `true`   | Emit MaveDB-ready CSVs per condition under `results/{exp}/deposit/mavedb/`.                                                                                                             |
| `run_cosmos`          | `false`  | Run cosmos multi-phenotype decomposition. Needs exactly two conditions with `phenotype` slots 1/2 in the experiment CSV. **Slow** (a model per position). Tune via the `cosmos:` block. |
| `run_qc`              | `true`   | Run FastQC + MultiQC.                                                                                                                                                                   |
| `noprocess`           | `false`  | Skip designed-variant filtering.                                                                                                                                                        |
| `remove_zeros`        | `false`  | Drop zero/unobserved-count variants before scoring.                                                                                                                                     |
| `regenerate_variants` | `false`  | Rebuild `variants_file` from `oligo_file`.                                                                                                                                              |
| `baseline_condition`  | —        | Condition used as the untreated/input baseline (library-quality QC).                                                                                                                    |
| `max_deletion_length` | `0`      | Max designed in-frame deletion (codons); longer insdel-resolved deletions rejected. `0` disables.                                                                                       |

### Optional `cosmos:` block

When `run_cosmos: true`, tune with `include_type` / `exclude_type` /
`x_gmm_n_components` / `min_num_variants_per_group`.

## Mapping and read processing

| Key                   | Default | Effect                                                                                                                                                                            |
| --------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `aligner`             | `bbmap` | `bbmap` (richer QC) or `minimap2` (much faster, slightly less QC).                                                                                                                |
| `kmers`               | `15`    | BBMap k-mer length.                                                                                                                                                               |
| `sam`                 | `"1.3"` | SAM version for BBMap (required for GATK compatibility).                                                                                                                          |
| `min_q`               | `30`    | GATK minimum base quality.                                                                                                                                                        |
| `min_variant_obs`     | `3`     | Minimum number of observations for GATK AnalyzeSaturationMutagenesis to include a variant.                                                                                        |
| `bbtools_compression` | `pigz`  | BBTools (de)compression: `pigz` (parallel), `bgzip`, or `none` (gzip; use if pigz absent or BBTools hangs). Legacy `bbtools_use_bgzip: true/false` still accepted with a warning. |
| `adapters`            | —       | Adapter FASTA for BBDuk trimming; `example.yaml` defaults to the bundled BBTools adapters.                                                                                                                                           |
| `contaminants`        | —       | Contaminant FASTAs for BBDuk removal; `example.yaml` defaults to PhiX + sequencing artifacts.                                                                                                                                       |

## Resources (memory)

Per-tool allocations in **MB**, one per heavy rule; each becomes that rule's
`resources: mem_mb`, which a cluster scheduler turns into `sbatch --mem` (Java
rules derive `-Xmx` a fixed headroom below). Defaults are sane — override only
to fit a specific machine/cluster.

| Key                  | Default (MB) | Rule                              |
| -------------------- | ------------ | --------------------------------- |
| `mem_bbduk`          | 2000         | trim/clean                        |
| `mem_bbmerge`        | 2000         | merge/correct                     |
| `mem_bbmap`          | 12000        | bbmap map + index (heaviest)      |
| `mem_minimap2`       | 1000         | minimap2 map + index              |
| `mem_gatk`           | 6000         | GATK AnalyzeSaturationMutagenesis |
| `mem_process_sample` | 2000         | per-sample variant processing     |
| `mem_fastqc`         | 1024         | FastQC                            |
| `mem_rosace`         | 16000        | Rosace scoring                    |
| `mem_rosace_aa`      | 16000        | rosace-aa scoring                 |
| `mem_lilace`         | 16000        | Lilace scoring                    |
| `mem_cosmos`         | 4000         | cosmos (run_cosmos)               |

`mem` (GB) is a legacy single BBTools-heap knob, superseded by the `mem_*`
allocations but retained for out-of-tree use.

## Local tool overrides

For platforms where the conda/renv environments don't resolve (e.g. ARM Macs),
point a rule at a locally-installed tool instead:

`samtools_local`, `rosace_local`, `lilace_local`, `rosace_aa_local` — all
`false` by default. See the main README "Troubleshooting" section.
