"""
Queyr sequence process
"""
from Bio import SeqIO
import PseudoBloomFilter

def reverse_com(string):
    rev_com = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    result = ""
    for i in string:
        result = rev_com[i] + result

    return result

class QuerySeq(object):
    """

    """
    def __init__(self, seq_record):
        self.seq = seq_record.seq
        self.id = seq_record.id
        self.length = len(self.seq)

    def check_kmer(self, bloom_filter):
        self.fw_hits = []
        self.rc_hits = []
        self.k = bloom_filter.get_k()
        num_bins = bloom_filter.get_num_bins()

        for i in range(self.length - self.k + 1):
            kmer = str(self.seq[i:i + self.k])
            reverse_kmer = reverse_com(kmer)
            for j in range(int(num_bins)):
                if bloom_filter.check_bin(kmer, j):
                    self.fw_hits.append((i, j, kmer))

                if bloom_filter.check_bin(reverse_kmer, j):
                    self.rc_hits.append((i, j, kmer))

if __name__ == '__main__':
    record1 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair2_query.fasta", "fasta")
    record2 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair2_target.fasta", "fasta")
    test_filter = PseudoBloomFilter.PseudoBloomFilter(record2, 9, 0, 0.75)
    test_filter.generate_filter()
    test_query = QuerySeq(record1)
    test_query.check_kmer(test_filter)
    print(test_query.fw_hits)

