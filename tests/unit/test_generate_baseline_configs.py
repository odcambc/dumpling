import pytest
import yaml

from workflow.rules.scripts.generate_baseline_configs import _run


def test_generate_baseline_configs_writes_expected_output(tmp_path, mock_snakemake):
    output_file = tmp_path / "baseline_config.yaml"
    multiqc_config = tmp_path / "multiqc_config.yaml"
    multiqc_config.write_text('extra_fn_clean_exts:\n  - ".existing"\n')

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

    result = yaml.safe_load(output_file.read_text())
    assert result["multiqc_dumpling"]["orf"] == "1-12"
    assert result["multiqc_dumpling"]["variants_file"] == (
        "config/designed_variants/example.csv"
    )
    assert ".existing" in result["extra_fn_clean_exts"]
    assert ".refCoverage" in result["extra_fn_clean_exts"]


def test_generate_baseline_configs_preserves_other_keys(tmp_path, mock_snakemake):
    """Unrelated top-level keys (custom_data, sp, etc.) must survive the merge."""
    output_file = tmp_path / "baseline_config.yaml"
    multiqc_config = tmp_path / "multiqc_config.yaml"
    multiqc_config.write_text(
        "extra_fn_clean_exts:\n"
        '  - "_trim"\n'
        "custom_data:\n"
        "  some_section:\n"
        '    file_format: "tsv"\n'
        "sp:\n"
        "  some_section:\n"
        '    fn: "*.tsv"\n'
    )

    snakemake = mock_snakemake(
        config={"orf": "1-12", "variants_file": "x.csv"},
        input={"multiqc_config": str(multiqc_config)},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    result = yaml.safe_load(output_file.read_text())
    assert result["custom_data"]["some_section"]["file_format"] == "tsv"
    assert result["sp"]["some_section"]["fn"] == "*.tsv"
    assert "_trim" in result["extra_fn_clean_exts"]
    assert ".refCoverage" in result["extra_fn_clean_exts"]


@pytest.mark.parametrize(
    "multiqc_yaml",
    [
        # Inline comment after the key — old line-equality match would miss this.
        'extra_fn_clean_exts:  # additional suffixes\n  - "_trim"\n',
        # CRLF line endings — old line-equality match would miss this.
        'extra_fn_clean_exts:\r\n  - "_trim"\r\n',
        # Trailing whitespace on the key line — old line-equality match would miss this.
        'extra_fn_clean_exts:   \n  - "_trim"\n',
    ],
    ids=["inline_comment", "crlf_line_endings", "trailing_whitespace"],
)
def test_generate_baseline_configs_handles_yaml_formatting_variants(
    tmp_path, mock_snakemake, multiqc_yaml
):
    """`.refCoverage` must be merged regardless of cosmetic YAML formatting."""
    output_file = tmp_path / "baseline_config.yaml"
    multiqc_config = tmp_path / "multiqc_config.yaml"
    multiqc_config.write_text(multiqc_yaml)

    snakemake = mock_snakemake(
        config={"orf": "1-12", "variants_file": "x.csv"},
        input={"multiqc_config": str(multiqc_config)},
        output=[str(output_file)],
        log=["/dev/null"],
    )

    _run(snakemake)

    result = yaml.safe_load(output_file.read_text())
    assert "_trim" in result["extra_fn_clean_exts"]
    assert ".refCoverage" in result["extra_fn_clean_exts"]


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
