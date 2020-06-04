import math
from pathlib import Path
from typing import List, Dict, Mapping, Union

import attr

from viral_verify.naive_bayes.constants import Classification
from viral_verify.naive_bayes.io import parse_naive_bayes_classifier_table, NaiveBayesClassifierFreqs


@attr.s
class NaiveBayesClassification:
    """Naive Bayes classification result"""
    contig_name: str = attr.ib()
    classification: str = attr.ib()
    log_viral_prob: float = attr.ib()
    log_plasmid_prob: float = attr.ib()
    log_chrom_prob: float = attr.ib()
    log_plasmid_or_chrom_prob: float = attr.ib()
    log_viral_minus_plasmid_or_chrom_prob: float = attr.ib()
    log_plasmid_minus_chrom_prob: float = attr.ib()
    uncertainty_threshold: float = attr.ib()

    @classmethod
    def from_contig_domains(cls,
                            contig: str,
                            domains: List[str],
                            classifier_table: Dict[str, NaiveBayesClassifierFreqs],
                            uncertainty_threshold: float = 3.0):
        log_chrom_prob = 0.0
        log_plasmid_prob = 0.0
        log_viral_prob = 0.0
        log_plasmid_or_chrom_prob = 0.0
        for domain in domains:
            try:
                freqs: NaiveBayesClassifierFreqs = classifier_table[domain]
                log_plasmid_prob += math.log(freqs.plasmid_freq)
                log_chrom_prob += math.log(freqs.chrom_freq)
                log_viral_prob += math.log(freqs.viral_freq)
                log_plasmid_or_chrom_prob += math.log(freqs.plasmid_or_chrom_freq)
            except KeyError:
                continue
        log_viral_minus_plasmid_or_chrom_prob = log_viral_prob - log_plasmid_or_chrom_prob
        log_plasmid_minus_chrom_prob = log_plasmid_prob - log_chrom_prob
        if log_viral_minus_plasmid_or_chrom_prob > uncertainty_threshold:
            classification = Classification.VIRUS
        elif log_viral_minus_plasmid_or_chrom_prob > (-1) * uncertainty_threshold:
            if len(domains) > 2:
                classification = Classification.UNCERTAIN_VIRAL_OR_BACTERIAL
            else:
                classification = Classification.UNCERTAIN_TOO_SHORT
        elif log_plasmid_minus_chrom_prob > uncertainty_threshold:
            classification = Classification.PLASMID
        elif (log_chrom_prob - log_plasmid_prob) > uncertainty_threshold:
            classification = Classification.CHROMOSOME
        else:
            classification = Classification.UNCERTAIN_PLASMID_OR_CHROMOSOMAL
        return cls(contig_name=contig,
                   classification=classification,
                   log_viral_prob=log_viral_prob,
                   log_plasmid_prob=log_plasmid_prob,
                   log_chrom_prob=log_chrom_prob,
                   log_plasmid_or_chrom_prob=log_plasmid_or_chrom_prob,
                   log_plasmid_minus_chrom_prob=log_plasmid_minus_chrom_prob,
                   uncertainty_threshold=uncertainty_threshold,
                   log_viral_minus_plasmid_or_chrom_prob=log_viral_minus_plasmid_or_chrom_prob)


def naive_bayes_classification(contig_domains: Mapping[str, List[str]],
                               classifier_table_path: Union[str, Path],
                               uncertainty_threshold: float = 3.0) -> Dict[str, NaiveBayesClassification]:
    """Naive Bayes classification of whether a contig is viral, plasmid or chromosomal.

    1. Calculate the probability a contig is viral, plasmid or chromosomal given the protein domains found in the contig
    using hmmsearch of Prodigal gene predictions against an HMM profile DB like Pfam.
    2. Given the frequencies that certain protein domains appear in viral, plasmid or chromosomal sequences, use Naive
    Bayesian method to classify sequences.
    """
    classifier_table: Dict[str, NaiveBayesClassifierFreqs] = parse_naive_bayes_classifier_table(classifier_table_path)
    out: Dict[str, NaiveBayesClassification] = {}
    for contig, domains in contig_domains.items():
        out[contig] = NaiveBayesClassification.from_contig_domains(contig=contig,
                                                                   domains=domains,
                                                                   classifier_table=classifier_table,
                                                                   uncertainty_threshold=uncertainty_threshold)
    return out
