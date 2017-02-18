from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


class Kmer(SeqRecord):
    def __init__(self, kmer_str, freq_dict=None):
        SeqRecord.__init__(self, Seq(kmer_str), id=kmer_str, description="")
        if freq_dict:
            self.freq = freq_dict.get(kmer_str, 1)
        else:
            self.freq = 1

    def check_threshold(self, threshold):
        if self.freq >= threshold:
            return True
        else:
            return False


    def to_fasta(self, handle, type):
        SeqIO.write(self, handle, type)

def main():
    return 0


if __name__ == '__main__':
    main()
