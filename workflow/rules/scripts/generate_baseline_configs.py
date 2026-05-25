import logging

import yaml

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

    with open(multiqc_config_file, "r") as f:
        cfg = yaml.safe_load(f) or {}

    extras = cfg.setdefault("extra_fn_clean_exts", [])
    if ".refCoverage" not in extras:
        extras.append(".refCoverage")

    output = {
        "multiqc_dumpling": {"orf": orf, "variants_file": variants_file},
    }
    output.update(cfg)

    logging.debug("Writing to file: %s", output_file)
    with open(output_file, "w") as f:
        yaml.safe_dump(output, f, sort_keys=False, default_flow_style=False)

    logging.debug("Finished writing to file: %s", output_file)


def main():
    run_script(snakemake, _run)


if __name__ == "__main__":
    main()
