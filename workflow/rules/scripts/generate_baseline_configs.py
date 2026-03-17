import logging

from snakemake.script import snakemake
from script_utils import run_script

"""Generate a MultiQC config file for the baseline samples. This will write the following
variables:
orf
extra_fn_clean_exts
variants_file
"""


def _run(snakemake):
    output_file = snakemake.output[0]
    multiqc_config_file = snakemake.input["multiqc_config"]

    orf = snakemake.config["orf"]
    variants_file = snakemake.config["variants_file"]

    extra_fn_clean_exts = [".refCoverage"]

    with open(output_file, "w+") as f:
        logging.debug("Writing to file: %s", output_file)
        f.write("multiqc_dumpling: " + "\n")
        f.write('  orf: "' + orf + '"\n')
        f.write('  variants_file: "' + variants_file + '"\n')
        f.write("\n")

        # Merge the shared MultiQC config into the generated baseline-specific config.
        with open(multiqc_config_file, "r") as multiqc_config:
            for line in multiqc_config:
                f.write(line)
                if line == "extra_fn_clean_exts:\n":
                    for ext in extra_fn_clean_exts:
                        f.write('  - "' + ext + '"\n')

    logging.debug("Finished writing to file: %s", output_file)


def main():
    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
