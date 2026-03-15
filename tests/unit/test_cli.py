from click.testing import CliRunner
from scanpy_workflow.cli import cli


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "scanpy-workflow" in result.output.lower() or "Usage" in result.output


def test_qc_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["qc", "--help"])
    assert result.exit_code == 0


def test_filter_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["filter", "--help"])
    assert result.exit_code == 0


def test_normalize_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["normalize", "--help"])
    assert result.exit_code == 0


def test_cluster_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["cluster", "--help"])
    assert result.exit_code == 0


def test_de_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["de", "--help"])
    assert result.exit_code == 0
