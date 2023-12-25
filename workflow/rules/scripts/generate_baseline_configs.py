import logging

"""Generate a MultiQC config file for the baseline samples. This will write the following
variables:
orf
extra_fn_clean_exts
variants_file
"""


output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)


orf = snakemake.config["orf"]
variants_file = snakemake.config["variants_file"]

extra_fn_clean_exts = [".refCoverage"]


with open(output_file, "w+") as f:
    logging.debug("Writing to file: %s", output_file)

    f.write('orf: "' + orf + '"\n')
    f.write('variants_file: "' + variants_file + '"\n')
    f.write("\n")

    # Now add the multiqc_config.yaml contents to the end of the file
    # We need to add additional extra_fn_clean_exts members to the list here
    with open("config/multiqc_config.yaml", "r") as multiqc_config:
        line = multiqc_config.readline()
        f.write(line)
        if line == "extra_fn_clean_exts:\n":
            for ext in extra_fn_clean_exts:
                f.write('- "' + ext + '"\n')

logging.debug("Finished writing to file: %s", output_file)
