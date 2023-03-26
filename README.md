# GATK based snakemake pipeline for deep mutational scanning experiments

This repository contains the snakemake-based workflow for implementing
deep mutational scanning experiments used in the Fraser and Coyote-Maestes labs.

Briefly, this conducts initial QC and mapping using bbtools, followed by the
AnalyzeSaturationMutagenesis GATK module to call variants in each replicate.