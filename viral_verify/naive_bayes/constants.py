from pkg_resources import resource_filename


class Classification:
    """Classification prediction constant values"""
    VIRUS = 'Virus'
    PLASMID = 'Plasmid'
    CHROMOSOME = 'Chromosome'
    UNCERTAIN_VIRAL_OR_BACTERIAL = 'Uncertain - viral or bacterial'
    UNCERTAIN_PLASMID_OR_CHROMOSOMAL = 'Uncertain - plasmid or chromosomal'
    UNCERTAIN_TOO_SHORT = 'Uncertain - too short'


CLASSIFIER_TABLE = resource_filename('viral_verify', 'data/classifier_table.txt')
DEFAULT_UNCERTAINTY_THRESHOLD = 3.0
