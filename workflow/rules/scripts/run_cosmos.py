"""Run cosmos (pimentellab/cosmos) on a dumpling-exported cosmos input CSV.

cosmos is a position-resolution causal model: given per-variant effects for an
*upstream* and a *downstream* phenotype (the `beta_hat_1`/`beta_hat_2` columns
that format_cosmos.py produces), it decomposes the downstream effect at each
position into the part mediated through the upstream phenotype and the part
that is direct. This is the "fully-fledged" cosmos stage: format_cosmos.py
prepares the wide input table, and this script *runs* cosmos on it and writes
the per-position decomposition.

The run flow (pinned against cosmos 1.0.2's vignette + API):

    data  = DMSData(df, phenotypes, include_type=..., exclude_type=...,
                    min_num_variants_per_group=...)
    prior = PriorFactory(data, x_name="beta_hat_1", y_name="beta_hat_2",
                         x_se_name="se_hat_1", x_gmm_n_components=...)
    model = ModelBuilder(prior, work_dir)
    for i in model.all_group_new_index:          # one position per iteration
        model.run_cosmos(i, no_s_hat=True, suppress_pareto_warning=True)
    summary = ModelAnalyzer(model, analysis_dir, has_position=True).summary(rank=1)

PERFORMANCE: run_cosmos is ~tens of seconds PER POSITION (a real optimization,
not a quick transform), and this script runs them serially. A large library
(hundreds of positions) takes hours. Parallelizing across positions
(scatter-gather) is a known future optimization — see tasks/tasks.md — but is
out of scope here; this is the straightforward single-process implementation.

First-cut scope: exactly two phenotypes (upstream -> downstream), i.e. two
conditions assigned phenotype slots 1 and 2. cosmos's x/y are the slot-1 and
slot-2 columns (beta_hat_1 -> beta_hat_2).
"""

import argparse
import tempfile

import pandas as pd

# Default cosmos type partitions. cosmos drops any variant whose `type` is in
# neither list, so this keeps missense and excludes synonymous + the
# length-resolved indel labels derive_cosmos_type emits (matching cosmos's own
# vignette exclude list).
DEFAULT_INCLUDE_TYPE = ["missense"]
DEFAULT_EXCLUDE_TYPE = [
    "synonymous",
    "nonsense",
    "insertion1",
    "insertion2",
    "insertion3",
    "deletion1",
    "deletion2",
    "deletion3",
]


def run_cosmos(
    input_csv,
    phenotypes,
    output_csv,
    include_type=None,
    exclude_type=None,
    x_gmm_n_components=2,
    min_num_variants_per_group=10,
):
    """Fit cosmos on a format_cosmos CSV and write the per-position summary.

    phenotypes : the two phenotype labels (slot order). Used by cosmos for its
                 phenotype count + plot labels; the causal x/y assignment is
                 fixed to beta_hat_1 -> beta_hat_2 by column position.
    """
    # Imported here (not at module top) so the module is importable without
    # cosmos installed; the cosmos conda env / the `cosmos` test extra provide it.
    from cosmos import DMSData, ModelAnalyzer, ModelBuilder, PriorFactory

    if len(phenotypes) != 2:
        raise ValueError(
            f"cosmos run requires exactly 2 phenotypes (upstream -> downstream); "
            f"got {len(phenotypes)}: {phenotypes}. Assign exactly two conditions "
            "phenotype slots 1 and 2 in the experiment CSV."
        )

    df = pd.read_csv(input_csv)
    data = DMSData(
        df,
        list(phenotypes),
        include_type=include_type or DEFAULT_INCLUDE_TYPE,
        exclude_type=exclude_type or DEFAULT_EXCLUDE_TYPE,
        min_num_variants_per_group=min_num_variants_per_group,
    )
    prior = PriorFactory(
        data,
        x_name="beta_hat_1",
        y_name="beta_hat_2",
        x_se_name="se_hat_1",
        x_gmm_n_components=x_gmm_n_components,
    )

    # cosmos writes per-position fit artifacts to disk; we only keep the final
    # summary, so point it at scratch dirs and emit just the summary CSV.
    with (
        tempfile.TemporaryDirectory() as work_dir,
        tempfile.TemporaryDirectory() as analysis_dir,
    ):
        model = ModelBuilder(prior, work_dir + "/")
        for group_idx in model.all_group_new_index:
            model.run_cosmos(group_idx, no_s_hat=True, suppress_pareto_warning=True)
        analyzer = ModelAnalyzer(model, analysis_dir + "/", has_position=True)
        summary = analyzer.summary(rank=1, save=False)

    summary.to_csv(output_csv, index=False)
    return summary


def main():
    # Snakemake injects `snakemake` as a module global; check globals() (same
    # gotcha documented in format_cosmos.main()).
    if "snakemake" in globals():
        cfg = snakemake.config.get("cosmos", {})
        run_cosmos(
            input_csv=snakemake.input.cosmos,
            phenotypes=list(snakemake.params.phenotypes),
            output_csv=snakemake.output.results,
            include_type=cfg.get("include_type"),
            exclude_type=cfg.get("exclude_type"),
            x_gmm_n_components=cfg.get("x_gmm_n_components", 2),
            min_num_variants_per_group=cfg.get("min_num_variants_per_group", 10),
        )
        return

    parser = argparse.ArgumentParser(description="Run cosmos on a cosmos input CSV.")
    parser.add_argument("--input", required=True, help="format_cosmos output CSV")
    parser.add_argument("--output", required=True, help="cosmos results CSV path")
    parser.add_argument(
        "--phenotypes",
        nargs=2,
        required=True,
        metavar=("UPSTREAM", "DOWNSTREAM"),
        help="the two phenotype labels, slot order (beta_hat_1 then beta_hat_2)",
    )
    parser.add_argument("--gmm-components", type=int, default=2)
    parser.add_argument("--min-variants-per-group", type=int, default=10)
    args = parser.parse_args()

    run_cosmos(
        input_csv=args.input,
        phenotypes=args.phenotypes,
        output_csv=args.output,
        x_gmm_n_components=args.gmm_components,
        min_num_variants_per_group=args.min_variants_per_group,
    )


if __name__ == "__main__":
    main()
