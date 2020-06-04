import attr
from Bio.SeqRecord import SeqRecord


@attr.s
class Contig:
    seq_rec: SeqRecord = attr.ib()
    is_circular: bool = attr.ib(default=False)
    n_matching_ends: int = attr.ib(default=0)
    min_length: int = attr.ib(default=500)
    kmax: int = attr.ib(default=200)
    kmin: int = attr.ib(default=50)

    @property
    def seq_len(self) -> int:
        return len(self.seq_rec.seq)

    def circular_seq(self) -> str:
        seq = str(self.seq_rec.seq)
        return seq if not self.is_circular else seq + seq[self.n_matching_ends:]

    def circular_seq_fasta(self) -> str:
        seq = self.circular_seq()
        return f'>{self.seq_rec.description}{" circular" if self.is_circular else ""}\n{seq}\n'

    @classmethod
    def from_seq_record(cls, rec: SeqRecord, min_length=500, kmax=200, kmin=50):
        seq_len = len(rec.seq)
        seq = str(rec.seq)
        is_circular, k = Contig.find_matching_at_ends(seq, kmax, kmin) if seq_len < min_length else False, 0
        return cls(seq_rec=rec,
                   is_circular=is_circular,
                   n_matching_ends=k,
                   min_length=min_length,
                   kmax=kmax,
                   kmin=kmin)

    @staticmethod
    def find_matching_at_ends(seq: str, kmax: int = 200, kmin: int = 50) -> (bool, int):
        """Find number of matching characters at the ends of a sequence

        Starting at a length of `kmax` (200) until `kmin` (50), see if there are up to `kmax` matching characters at
        the ends of `seq`.
        """
        k = 0
        is_circular = False
        seq_len = len(seq)
        for k in range(kmax, kmin, -1):
            if k >= seq_len:
                continue
            start: str = seq[:k]
            end: str = seq[-k:]
            if start == end:
                is_circular = True
                break
        return is_circular, k
