# -*- coding: utf-8 -*-
"""Constants for hmmsearch module"""
import re

REGEX_PRODIGAL_GENE_NUMBER = re.compile(r'_\d+$')
"""Regular expression to match Prodigal gene number appended to contig sequence ID"""

HMMSEARCH_DOMTBLOUT_COLUMNS = ['target_name',
                               'target_accession',
                               'target_length',
                               'query_name',
                               'query_accession',
                               'query_length',
                               'full_sequence_evalue',
                               'full_sequence_score',
                               'full_sequence_bias',
                               'domain_number',
                               'domain_of',
                               'domain_c_evalue',
                               'domain_i_evalue',
                               'domain_score',
                               'domain_bias',
                               'hmm_coord_from',
                               'hmm_coord_to',
                               'ali_coord_from',
                               'ali_coord_to',
                               'env_coord_from',
                               'env_coord_to',
                               'env_coord_acc', ]
"""HMMer3 hmmsearch domtblout output table columns except last column containing free text of sequence description"""
