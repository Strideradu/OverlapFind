from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class Kmer(SeqRecord):
    def __init__(self, kmer_str, freq_dict=None):
        SeqRecord.__init__(Seq(kmer_str), id=kmer_str, description="")
        if freq_dict:
            self.freq = freq_dict[kmer_str]
        else:
            self.freq = 1

    def check_threshold(self, threshold):
        if self.freq >= threshold:
            return True
        else:
            return False


def main():
    return 0


if __name__ == '__main__':
    main()
