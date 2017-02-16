from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class kmer(SeqRecord):
    def __init__(self, kmer_str, freq_dict = None):
        SeqRecord.__init__(Seq(kmer_str), id=kmer_str, description="")
        if freq_dict:
            self.freq = freq_dict[kmer_seq]
        else:
            self.freq = 1

def main():
    return 0

if __name__ == '__main__':
    main()
