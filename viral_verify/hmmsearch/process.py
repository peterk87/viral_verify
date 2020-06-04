import logging
import subprocess as sp
from pathlib import Path
from typing import Union, IO

logger = logging.getLogger(__name__)


def run_hmmsearch(hmm_db: Union[str, Path, IO],
                  input_fasta: Union[str, Path, IO],
                  raw_output: Union[str, Path, IO],
                  tblout: Union[str, Path, IO],
                  threads: int = 1) -> None:
    """Run HMMer3 hmmsearch with a protein sequence FASTA against an HMM profile DB like Pfam.

    Parameters
    ----------
    hmm_db : Union[str, Path, IO]
        HMM profile DB path
    input_fasta : Union[str, Path, IO]
        Protein FASTA path
    raw_output : Union[str, Path, IO]
        Raw hmmsearch output path
    tblout : Union[str, Path, IO]
        Tabular space-delimited ``--domtblout`` protein domain output path
    threads : int
        Number of threads to run hmmsearch with
    """
    cmd_list = ['hmmsearch', '--noali', '--cut_nc',
                '-o', str(raw_output), '--domtblout', str(tblout),
                '--cpu', str(threads),
                str(hmm_db), str(input_fasta)]
    cmd = ' '.join(cmd_list)
    logger.info(f'Running hmmsearch command: {cmd}')
    sp.run(cmd_list,
           stdout=sp.PIPE,
           stderr=sp.PIPE,
           check=True)
    logger.info(f'Ran hmmsearch with tabular output at "{tblout}" and raw output at "{raw_output}"')
