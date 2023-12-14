import pandas as pd
import logging

"""Generate a list of files specific to baseline conditionsto be used as input for MultiQC.
This is just a list of files that will be passed using --file-list

The baseline QC reads the following files:
- fastqc
- bbduk
- bbmap
- gatk analyzesaturationmutagenesis
- count processing script

Essentially, it does not continue to Enrich2, since we don't run that on baseline samples.
"""


output_file = snakemake.output[0]
log_file = snakemake.log[0]
experiment_name = snakemake.config["experiment"]

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)

logging.debug("Writing to file: %s", output_file)

baseline_condition = snakemake.config["baseline_condition"]

experiments = pd.read_csv(snakemake.config["experiment_file"], header=0).dropna(how = 'all').set_index(
    "sample", drop=False, verify_integrity=True
)

baseline_samples = experiments.loc[
    (experiments["condition"] == baseline_condition)
]["sample"]

with open(output_file, "w+") as f:
    logging.debug("Writing to file: %s", output_file)

    for sample in baseline_samples.index:
        # Write fastqc output (Assuming paired-end reads)
        f.write("./stats/" + experiment_name + "/fastqc/" + sample + "R1_001_fastqc.html\n")
        f.write("./stats/" + experiment_name + "/fastqc/" + sample + "R2_001_fastqc.html\n")
        # Write bbduk output
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.qhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.bhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.gchist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.aqhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.lhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim.stats.txt\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_trim_contam.stats.txt\n")
        # Write bbmerge output
        f.write("./stats/" + experiment_name + "/" + sample + "_merge.ihist\n")
        # Write bbmap output
        f.write("./stats/" + experiment_name + "/" + sample + "_map.covstats\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.basecov\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.bincov\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.ehist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.indelhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.mhist\n")
        f.write("./stats/" + experiment_name + "/" + sample + "_map.idhist\n")
        # Write GATK ASM output
        f.write("./results/" + experiment_name + "/gatk/" + sample + ".coverageLengthCounts\n")
        f.write("./results/" + experiment_name + "/gatk/" + sample + ".readCounts\n")
        f.write("./results/" + experiment_name + "/gatk/" + sample + ".refCoverage\n")
        # Write count processing script output
        f.write("./results/" + experiment_name + "/processed_counts/counts/" + sample + ".csv\n")

logging.debug("Done")