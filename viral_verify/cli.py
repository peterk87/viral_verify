"""Console script for viral_verify."""
import logging
import multiprocessing
import sys
from pathlib import Path
from typing import Optional, Dict

import click

from viral_verify.contig import Contig
from viral_verify.hmmsearch import top_hmm_results, run_hmmsearch
from viral_verify.io import parse_contigs, write_circular_contigs_fasta, filter_predicted_genes, \
    output_classified_contigs, output_results_table, hmm_names_to_desc
from viral_verify.log import init_logging
from viral_verify.naive_bayes import naive_bayes_classification, DEFAULT_UNCERTAINTY_THRESHOLD, CLASSIFIER_TABLE
from viral_verify.prodigal import prodigal_meta

logger = logging.getLogger(__name__)


@click.command()
@click.option('-i', '--input-fasta', type=click.Path(exists=True),
              required=True, help='Input fasta file')
@click.option('-o', '--outdir', type=click.Path(exists=False, writable=True),
              required=True, help='Output directory')
@click.option('-H', '--hmm-db', type=click.Path(exists=True),
              required=True, help='Path to Pfam-A HMM database')
@click.option('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
              help=f'Number of threads (default={multiprocessing.cpu_count()})')
@click.option('-p', '--output-plasmids-separately', is_flag=True,
              help='Output predicted plasmids separately?')
@click.option('--prefix', default=None, help='Output file prefix (default: None)')
@click.option('--uncertainty-threshold', type=float, default=DEFAULT_UNCERTAINTY_THRESHOLD,
              help=f'Uncertainty threshold (Natural log probability) (default={DEFAULT_UNCERTAINTY_THRESHOLD})')
@click.option('--naive-bayes-classifier-table', type=click.Path(exists=True), default=CLASSIFIER_TABLE,
              help=f'Table of protein domain frequencies to use for Naive Bayes classification '
                   f'(default="{CLASSIFIER_TABLE}")')
@click.option('-v', '--verbose', default=2, count=True, help='Logging verbosity')
@click.version_option()
def main(input_fasta: str,
         outdir: str,
         hmm_db: str,
         threads: int,
         output_plasmids_separately: bool,
         prefix: Optional[str],
         uncertainty_threshold: float,
         naive_bayes_classifier_table: str,
         verbose: int):
    """HMM and Naive Bayes classification of contig sequences as either viral, plasmid or chromosomal.

    Requires Prodigal for gene prediction and hmmsearch from HMMer3 for searching for Pfam HMM profiles.
    """
    init_logging(verbose)
    input_fasta_path = Path(input_fasta).resolve()
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True)
    if prefix:
        logger.info(f'Output file prefix="{prefix}"')
    else:
        prefix = input_fasta_path.stem
        logger.info(f'Output file prefix not specified. Using input FASTA filename as prefix ("{prefix}")')

    logger.info(f'Parsing contig sequences from "{input_fasta}" and determine if any could be circular')
    contig_infos: Dict[str, Contig] = parse_contigs(input_fasta)
    logger.info(f'Parsed {len(contig_infos)} contigs from "{input_fasta}"')
    input_fasta_circularized: Path = outdir_path / (prefix + "-circularized.fasta")
    logger.info(f'Writing circularized contig sequences to "{input_fasta_circularized}"')
    write_circular_contigs_fasta(contig_infos, input_fasta_circularized)
    logger.info(f'Running Prodigal gene prediction on "{input_fasta_circularized}"')
    proteins_fasta_path: Path = outdir_path / (prefix + '-proteins.fa')
    genes_fasta_path: Path = outdir_path / (prefix + '-genes.fa')
    prodigal_meta(input_fasta=input_fasta_circularized,
                  genes_fasta=genes_fasta_path,
                  proteins_fasta=proteins_fasta_path)
    logger.info(f'Prodigal protein sequence output at "{proteins_fasta_path}"')
    logger.info(f'Prodigal nucleotide sequence output at "{genes_fasta_path}"')

    filtered_proteins_path = outdir_path / (prefix + "-proteins-circularized.fa")
    logger.info(f'Filtering out genes predicted over the end of the expected end of each contig. '
                f'Output at "{filtered_proteins_path}"')
    filter_predicted_genes(proteins_fasta_path, filtered_proteins_path, contig_infos)
    logger.info(f'Parsing HMM names and descriptions from "{hmm_db}"')
    protein_name_to_desc = hmm_names_to_desc(hmm_db)
    logger.info(f'Parsed {len(protein_name_to_desc)} names and descriptions from "{hmm_db}"')
    logger.info(f'hmmsearch of "{filtered_proteins_path}" against HMM DB "{hmm_db}" with {threads} threads.')
    hmmsearch_raw_output = outdir_path / (prefix + '-hmmsearch.output')
    hmmsearch_tblout = outdir_path / (prefix + '-hmmsearch.domtblout')
    run_hmmsearch(hmm_db=hmm_db,
                  input_fasta=filtered_proteins_path,
                  raw_output=hmmsearch_raw_output,
                  tblout=hmmsearch_tblout,
                  threads=threads)
    logger.info(f'hmmsearch raw results output at "{hmmsearch_raw_output}"')
    logger.info(f'hmmsearch tabular output at "{hmmsearch_tblout}"')
    logger.info(f'Parsing hmmsearch tabular output "{hmmsearch_tblout}"')
    contig_domains, top_domains = top_hmm_results(hmmsearch_tblout)
    logger.info(f'Parsed {sum(1 for k, vs in top_domains.items() for v in vs)} protein domain results for '
                f'{len(top_domains)} contigs (out of {len(contig_infos)} total contigs) from hmmsearch tabular '
                f'output "{hmmsearch_tblout}"')
    logger.debug(f'Top domains={top_domains}')
    logger.debug(f'Contig domains={contig_domains}')
    logger.info(f'Naive Bayes classification of top predicted protein domains in each contig')
    contig_classifications = naive_bayes_classification(contig_domains=contig_domains,
                                                        classifier_table_path=naive_bayes_classifier_table,
                                                        uncertainty_threshold=uncertainty_threshold)
    results_csv_path = outdir_path / (prefix + '-results.csv')
    logger.info(f'Writing output results CSV to "{results_csv_path}"')
    output_results_table(results_csv_path=results_csv_path,
                         contigs=contig_infos,
                         contig_domains=contig_domains,
                         contig_classifications=contig_classifications,
                         protein_name_to_desc=protein_name_to_desc)
    output_classified_contigs(contig_classifications=contig_classifications,
                              contigs=contig_infos,
                              outdir=outdir_path,
                              output_plasmids_separately=output_plasmids_separately,
                              prefix=prefix)

    logger.info(f'Done! Results can be found in "{outdir_path}". '
                f'Classification results can be found at "{results_csv_path}"')


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
