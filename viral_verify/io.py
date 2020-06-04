import re
from pathlib import Path
from typing import Dict, Union, IO, List, Mapping, Iterator

import attr
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from viral_verify.contig import Contig
from viral_verify.naive_bayes import NaiveBayesClassification
from viral_verify.naive_bayes.constants import Classification
from viral_verify.prodigal import prodigal_gene_start


def write_circular_contigs_fasta(contig_infos: Dict[str, Contig],
                                 output_path: Union[str, Path, IO]) -> None:
    with open(output_path, 'w') as fout:
        for contig_info in contig_infos.values():
            fout.write(contig_info.circular_seq_fasta())


def parse_contigs(fasta_path: Union[str, Path, IO],
                  min_length: int = 500,
                  kmax: int = 200,
                  kmin: int = 50) -> Dict[str, Contig]:
    return {rec.id: Contig.from_seq_record(rec,
                                           min_length=min_length,
                                           kmax=kmax,
                                           kmin=kmin) for rec in SeqIO.parse(fasta_path, 'fasta')}


def output_classified_contigs(contig_classifications: Dict[str, NaiveBayesClassification],
                              contigs: Dict[str, Contig],
                              outdir: Path,
                              output_plasmids_separately: bool,
                              prefix: str) -> None:
    prediction_fasta_dir: Path = outdir / 'classified-fasta-output'
    prediction_fasta_dir.mkdir(parents=True, exist_ok=True)
    viral_recs = []
    plasmid_recs = []
    chromosome_recs = []
    viral_uncertain_recs = []
    plasmid_uncertain_recs = []
    unclassified_recs = []

    contig: str
    info: Contig
    for contig, info in contigs.items():
        if contig not in contig_classifications:
            continue
        nbc: NaiveBayesClassification = contig_classifications[contig]
        # output_fasta_entry = f'>{info.seq_rec.description}\n{str(info.seq_rec.seq)}\n'
        if nbc.classification == Classification.VIRUS:
            viral_recs.append(info.seq_rec)
        elif nbc.classification == Classification.CHROMOSOME:
            chromosome_recs.append(info.seq_rec)
        elif nbc.classification == Classification.PLASMID:
            if output_plasmids_separately:
                plasmid_recs.append(info.seq_rec)
            else:
                chromosome_recs.append(info.seq_rec)
        elif nbc.classification == Classification.UNCERTAIN_VIRAL_OR_BACTERIAL:
            viral_uncertain_recs.append(info.seq_rec)
        elif nbc.classification == Classification.UNCERTAIN_PLASMID_OR_CHROMOSOMAL:
            if output_plasmids_separately:
                plasmid_uncertain_recs.append(info.seq_rec)
            else:
                chromosome_recs.append(info.seq_rec)
        else:
            unclassified_recs.append(info.seq_rec)

    if len(viral_recs) > 0:
        SeqIO.write(viral_recs,
                    prediction_fasta_dir / (prefix + '-viral.fasta'),
                    'fasta')
    if len(plasmid_recs) > 0:
        SeqIO.write(plasmid_recs,
                    prediction_fasta_dir / (prefix + '-plasmid.fasta'),
                    'fasta')
    if len(chromosome_recs) > 0:
        SeqIO.write(chromosome_recs,
                    prediction_fasta_dir / (prefix + '-chromosome.fasta'),
                    'fasta')
    if len(viral_uncertain_recs) > 0:
        SeqIO.write(viral_uncertain_recs,
                    prediction_fasta_dir / (prefix + '-viral_uncertain.fasta'),
                    'fasta')
    if len(plasmid_uncertain_recs) > 0:
        SeqIO.write(plasmid_uncertain_recs,
                    prediction_fasta_dir / (prefix + '-plasmid_uncertain.fasta'),
                    'fasta')
    if len(unclassified_recs) > 0:
        SeqIO.write(unclassified_recs,
                    prediction_fasta_dir / (prefix + '-unclassified.fasta'),
                    'fasta')


def filter_predicted_genes(input_fasta: Union[str, Path, IO],
                           output_fasta: Union[str, Path, IO],
                           contig_len_circ: Dict[str, Contig]) -> None:
    filtered_recs: List[SeqRecord] = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        contig_name = re.sub(r'_\d+$', '', rec.id)
        gene_start = prodigal_gene_start(rec.description)
        contig_info = contig_len_circ[contig_name]
        if gene_start < contig_info.seq_len:
            filtered_recs.append(rec)
    SeqIO.write(filtered_recs, output_fasta, 'fasta')


def output_results_table(results_csv_path: Path,
                         contigs: Dict[str, Contig],
                         contig_domains: Mapping[str, List[str]],
                         contig_classifications: Dict[str, NaiveBayesClassification],
                         protein_name_to_desc: Dict[str, str]) -> None:
    results: List[Dict] = []
    for contig_name, contig in contigs.items():
        if contig_name in contig_classifications:
            result = attr.asdict(contig_classifications[contig_name])
            result['protein_domains'] = ';'.join((f'{x} [{protein_name_to_desc[x]}]'
                                                  for x in contig_domains[contig_name]))
            results.append(result)
        else:
            results.append(dict(contig_name=contig_name,
                                classification='Unclassified'))
    df_results = pd.DataFrame(results)
    df_results.to_csv(results_csv_path, index=False)


def parse_hmms(hmm_path: Union[str, Path]) -> Iterator[str]:
    """Simple HMM file parser from HMMer3 build

    Parameters
    ----------
    hmm_path
        Path to HMM file build by HMMer3

    Yields
    ------
    str
        HMM entry text
    """
    with open(hmm_path) as handle:
        entry = ''
        for line in handle:
            if line.startswith('//'):
                entry += line
                yield entry
                entry = ''
            elif line.rstrip() == '':
                continue
            else:
                entry += line
        if entry != '':
            yield entry


def hmm_name(entry_text: str) -> str:
    """Get the HMM NAME attribute value"""
    for line in entry_text.splitlines():
        if line.startswith('NAME'):
            return line.replace('NAME', '').strip()


def hmm_desc(entry_text: str) -> str:
    """Get the HMM DESC attribute value"""
    for line in entry_text.splitlines():
        if line.startswith('DESC'):
            return line.replace('DESC', '').strip()


def hmm_names_to_desc(hmm_path: Union[str, Path]) -> Dict[str, str]:
    """Parse an HMM file into a dictionary of protein names to descriptions"""
    return {hmm_name(x): hmm_desc(x) for x in parse_hmms(hmm_path)}
