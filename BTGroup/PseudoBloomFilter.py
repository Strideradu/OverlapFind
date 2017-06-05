"""
Process the ref sequence
"""
import ProbFunc
from Bio import SeqIO


class PseudoBloomFilter(object):
    def __init__(self, seq_record, k, L=0, match_rate=0.0):
        """

        :param seq: a biopython seq object
        :param k: size of k
        :param L: group hit distance threshold
        :param match_rate: sequence match rate
        """
        if L != 0:
            self.L = L
        else:
            if match_rate == 0:
                raise ValueError

            else:
                self.L = ProbFunc.statistical_bound_of_waiting_time(match_rate, k)

        self.seq = seq_record.seq
        self.id = seq_record.id
        self.length = len(self.seq)
        self.k = k

    def generate_filter(self):
        """
        generate bloom filter
        :return:
        """
        self.num_bins = self.length // self.L + 1

        bin = [set() for row in range(self.num_bins)]

        for i in range(self.length - self.k + 1):
            kmer = str(self.seq[i:i + self.k])
            bin_index = i // self.L
            bin[bin_index].add(kmer)

        self.bin = bin

    def check_bin(self, kmer, bin_index):
        return kmer in self.bin[bin_index]

    def get_num_bins(self):

        return self.num_bins

    def get_k(self):
        return self.k


if __name__ == '__main__':
    record1 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair2_query.fasta", "fasta")
    record2 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair2_target.fasta", "fasta")
    test_filter = PseudoBloomFilter(record1, 9, 0, 0.75)
    print test_filter.generate_filter()
