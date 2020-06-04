# -*- coding: utf-8 -*-
"""I/O functions for hmmsearch"""
import collections
from collections import OrderedDict
from pathlib import Path
from typing import Union, IO, List, Tuple, Mapping, Dict

import pandas as pd

from viral_verify.hmmsearch.constants import REGEX_PRODIGAL_GENE_NUMBER
from viral_verify.hmmsearch.result import HmmSearchResult


def parse_domtblout(domtblout: Union[str, Path, IO]) -> List[HmmSearchResult]:
    """Parse an HMMer3 hmmsearch domtblout table of protein domain predictions into a list of HmmSearchResult"""
    with open(domtblout) as fh:
        out = []
        for line in fh:
            if line.startswith('#'):
                continue
            out.append(HmmSearchResult.from_row_str(line))

    def sort_key(hmm_result: HmmSearchResult):
        return hmm_result.target_name, hmm_result.domain_score, hmm_result.ali_coord_from

    out.sort(key=sort_key, reverse=True)
    return out


def domtblout_to_dataframe(tblout: Union[str, Path, IO]) -> pd.DataFrame:
    """Parse an HMMer3 hmmsearch domtblout table of protein predictions into a Pandas DataFrame"""
    from viral_verify.hmmsearch.constants import HMMSEARCH_DOMTBLOUT_COLUMNS
    df = pd.read_table(tblout,
                       sep=r'\s+',
                       header=None,
                       comment='#',
                       usecols=range(0, 22),
                       names=HMMSEARCH_DOMTBLOUT_COLUMNS)
    df.sort_values(['target_name', 'domain_score', 'ali_coord_from'], ascending=False, inplace=True)
    return df


def top_hmm_results(tblout: Union[str, Path, IO]) -> Tuple[Mapping[str, List[str]], Dict[str, List[HmmSearchResult]]]:
    """Get the top hmmsearch predicted domains for each gene and contig

    Parameters
    ----------
    tblout
        Tabular space-delimited hmmsearch output path (``--domtblout``)

    Returns
    -------
    Tuple[Mapping[str, List[str]], Dict[str, List[HmmSearchResult]]]
        2 element tuple: dict of contig name to top predicted protein domains; dict of Prodigal gene name
        to list of HmmSearchResult
    """
    hmm_results: List[HmmSearchResult] = parse_domtblout(tblout)
    top_domains = top_domains_per_predicted_gene(hmm_results)
    contig_domains = top_domains_per_contig(top_domains)
    return contig_domains, top_domains


def top_domains_per_contig(top_domains: Dict[str, List[HmmSearchResult]]) -> Mapping[str, List[str]]:
    contig_domains: OrderedDict[str, List[str]] = collections.OrderedDict()
    results: List[HmmSearchResult]
    for target_name, results in top_domains.items():
        contig_name: str = REGEX_PRODIGAL_GENE_NUMBER.sub('', target_name)
        if contig_name not in contig_domains:
            contig_domains[contig_name] = [r.query_name for r in results]
        else:
            contig_domains[contig_name] += [r.query_name for r in results]
    return contig_domains


def top_domains_per_predicted_gene(hmm_results: List[HmmSearchResult]) -> Dict[str, List[HmmSearchResult]]:
    top_domains = {}
    for hmm_result in hmm_results:
        target_name: str = hmm_result.target_name
        if target_name not in top_domains:
            top_domains[target_name] = [hmm_result]
        else:
            if not hmm_result.overlaps(top_domains[target_name]):
                top_domains[target_name] += [hmm_result]
    return top_domains
