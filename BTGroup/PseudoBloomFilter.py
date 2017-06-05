"""
Process the ref sequence
"""
import ProbFunc

class PseudoBloomFilter(object):
    def __init__(self, seq_record, k, L=0, match_rate = 0):
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
        num_bins = self.length//self.L + 1

        bin = [set() for row in range(num_bins)]

        for i in range(self.length - self.k + 1):
            kmer = str(self.seq[i:i+self.k])
            bin_index = i//self.L
            bin[bin_index].add(kmer)
