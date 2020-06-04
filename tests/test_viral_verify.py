#!/usr/bin/env python

"""Tests for `viral_verify` package."""
import subprocess
from pathlib import Path

import pandas as pd
from click.testing import CliRunner
from pandas.testing import assert_frame_equal

from viral_verify import cli


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert 'viral_verify' in help_result.stdout

    result = runner.invoke(cli.main)
    assert result.exit_code == 2
    assert 'Missing option' in result.output

    test_fasta = Path('tests/data/test.fasta').resolve()
    hmm_db_gz = Path('tests/data/Pfam-A-filtered-for-tests.hmm.gz').resolve()
    expected_results_csv_path = Path('tests/data/expected-viral_verify-results.csv').resolve()
    with runner.isolated_filesystem():
        outdir = 'outdir'
        hmm_db = 'Pfam-A-filtered-for-tests.hmm'
        with open(hmm_db, 'w') as f:
            subprocess.call(['gunzip', '-c', hmm_db_gz.absolute()], stdout=f)

        result = runner.invoke(cli.main, ['-i', test_fasta.absolute(),
                                          '-o', outdir,
                                          '--hmm-db', hmm_db,
                                          '--prefix', 'test',
                                          '-p',
                                          '-vvv'])
        assert result.exit_code == 0
        outdir_path = Path(outdir)
        assert outdir_path.exists()
        assert (outdir_path / 'classified-fasta-output' / 'test-viral.fasta').exists()
        assert (outdir_path / 'classified-fasta-output' / 'test-unclassified.fasta').exists()
        assert (outdir_path / 'classified-fasta-output' / 'test-plasmid.fasta').exists()

        results_csv_path = outdir_path / 'test-results.csv'
        df_observed = pd.read_csv(results_csv_path)
        df_expected = pd.read_csv(expected_results_csv_path)
        assert_frame_equal(df_observed, df_expected)
