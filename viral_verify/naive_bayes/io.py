from typing import Dict

import attr


@attr.s
class NaiveBayesClassifierFreqs:
    """Naive Bayes classification table frequency values entry

    NaiveBayesClassifierFreqs contains the frequency values for whether a particular domain is found in viral, plasmid
    or chromosomal sequence."""
    name: str = attr.ib()
    plasmid_freq: float = attr.ib()
    chrom_freq: float = attr.ib()
    viral_freq: float = attr.ib()
    plasmid_or_chrom_freq: float = attr.ib()

    @classmethod
    def from_line(cls, line):
        name, _, _, _, plasmid_freq, chrom_freq, viral_freq, plasmid_or_chrom_freq = line.strip().split('\t')
        return cls(name=name,
                   plasmid_freq=float(plasmid_freq),
                   chrom_freq=float(chrom_freq),
                   viral_freq=float(viral_freq),
                   plasmid_or_chrom_freq=float(plasmid_or_chrom_freq))


def parse_naive_bayes_classifier_table(path) -> Dict[str, NaiveBayesClassifierFreqs]:
    out = {}
    with open(path) as f:
        for line in f:
            freqs = NaiveBayesClassifierFreqs.from_line(line)
            out[freqs.name] = freqs
    return out
