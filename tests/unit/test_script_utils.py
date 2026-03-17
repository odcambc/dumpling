import pytest

from workflow.rules.scripts.script_utils import run_script


def test_run_script_logs_and_reraises(mock_snakemake, mocker):
    snakemake = mock_snakemake(log=["/dev/null"])
    basic_config = mocker.patch("workflow.rules.scripts.script_utils.logging.basicConfig")
    log_exception = mocker.patch("workflow.rules.scripts.script_utils.logging.exception")

    def boom(_snakemake):
        raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        run_script(snakemake, boom)

    basic_config.assert_called_once()
    log_exception.assert_called_once_with("Script failed")
