from typing import List

import attr


@attr.s
class HmmSearchResult:
    """HMMer3 hmmsearch protein domain prediction result

    Only values necessary for Naive Bayes classification are used, i.e. domain
    score and alignment coordinates for checking overlap with other predicted
    domains."""
    target_name: str = attr.ib(init=True)
    """Target sequence name (e.g. predicted protein gene ID)"""
    query_name: str = attr.ib(init=True)
    """Query sequence name (e.g. protein domain name)"""
    domain_score: float = attr.ib(init=True, converter=float, validator=attr.validators.instance_of(float))
    """hmmsearch protein domain prediction score"""
    ali_coord_from: int = attr.ib(init=True, converter=int, validator=attr.validators.instance_of(int))
    """hmmsearch result alignment start index"""
    ali_coord_to: int = attr.ib(init=True, converter=int, validator=attr.validators.instance_of(int))
    """hmmsearch result alignment end index"""

    @classmethod
    def from_row_str(cls, row: str):
        """Parse HmmSearchResult from a non-comment (#) line of hmmsearch domtblout table output"""
        cells = row.strip().split()
        return cls(target_name=cells[0],
                   query_name=cells[3],
                   domain_score=float(cells[13]),
                   ali_coord_from=int(cells[17]),
                   ali_coord_to=int(cells[18]))

    def overlaps(self, target_results: List['HmmSearchResult']) -> bool:
        """Does this HmmSearchResult overlap other HmmSearchResults?"""
        return any((self.ali_coord_to > t.ali_coord_from and self.ali_coord_from < t.ali_coord_to
                    for t in target_results))
