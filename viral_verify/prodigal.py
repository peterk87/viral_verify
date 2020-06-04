import logging
import subprocess as sp
from pathlib import Path
from typing import Union, IO

logger = logging.getLogger(__name__)


def prodigal_meta(input_fasta: Union[str, Path, IO],
                  genes_fasta: Union[str, Path, IO],
                  proteins_fasta: Union[str, Path, IO]) -> None:
    """Run Prodigal gene prediction in metagenomic mode on an input FASTA format file

    Produces files containing nucleotide and protein sequences of Prodigal predicted genes
    """
    cmd_list = ['prodigal', '-p', 'meta', '-c',
                '-i', str(input_fasta),
                '-a', str(proteins_fasta),
                '-o', str(genes_fasta)]
    cmd = " ".join(cmd_list)
    logger.info(f'Running Prodigal gene prediction in metagenomic mode with command: {cmd}')
    sp.run(cmd_list,
           stdout=sp.PIPE,
           stderr=sp.PIPE,
           check=True)
    logger.info(f'Ran Prodigal gene prediction outputting protein sequences to '
                f'"{proteins_fasta}" and nucleotide sequences to  "{genes_fasta}".')


def prodigal_gene_start(rec_description: str) -> int:
    """Get a gene start index from a Prodigal FASTA header

    Examples
    --------
    Given the following Prodigal FASTA output header, parse the gene start index (i.e. 197)
    >>> prodigal_gene_start("k141_2229_1 # 197 # 379 # 1 # ID=4_1;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.437")
    197

    Parameters
    ----------
    rec_description : str
        SeqRecord description of Prodigal FASTA header

    Returns
    -------
    int
        Gene start index
    """
    return int(rec_description.split('#')[1].strip())
