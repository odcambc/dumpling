import pytest

from workflow.rules.scripts.generate_baseline_configs import _run


def test_generate_baseline_configs_writes_expected_output(tmp_path, mock_snakemake):
    output_file = tmp_path / "baseline_config.yaml"
    multiqc_config = tmp_path / "multiqc_config.yaml"
    multiqc_config.write_text("extra_fn_clean_exts:\n  - \".existing\"\n")

    snakemake = mock_snakemake(
        config={
            "orf": "1-12",
            "variants_file": "config/designed_variants/example.csv",
        },
        input={"multiqc_config": str(multiqc_config)},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    output = output_file.read_text()
    assert 'orf: "1-12"' in output
    assert 'variants_file: "config/designed_variants/example.csv"' in output
    assert '  - ".existing"' in output
    assert '  - ".refCoverage"' in output


def test_generate_baseline_configs_missing_multiqc_input(tmp_path, mock_snakemake):
    output_file = tmp_path / "baseline_config.yaml"
    missing_multiqc_config = tmp_path / "missing_multiqc_config.yaml"

    snakemake = mock_snakemake(
        config={
            "orf": "1-12",
            "variants_file": "config/designed_variants/example.csv",
        },
        input={"multiqc_config": str(missing_multiqc_config)},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    with pytest.raises(FileNotFoundError):
        _run(snakemake)
