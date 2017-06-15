"""
Queyr sequence process
"""
from Bio import SeqIO
import PseudoBloomFilter
from matplotlib import pyplot as plt

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
        self.L = bloom_filter.get_L()
        self.target_length = bloom_filter.get_length()
        num_bins = bloom_filter.get_num_bins()

        for i in range(self.length - self.k + 1):
            kmer = str(self.seq[i:i + self.k])
            reverse_kmer = reverse_com(kmer)
            for j in range(int(num_bins)):
                fw_pos = bloom_filter.check_bin(kmer, j)
                if fw_pos:
                    self.fw_hits.append((i, j, kmer, fw_pos))

                rc_pos = bloom_filter.check_bin(reverse_kmer, j)
                if rc_pos:
                    self.rc_hits.append((i, j, kmer, rc_pos))

        self.fw_hits.sort()
        self.rc_hits.sort()

    def fw_diag_group(self):
        best_cluster = None
        best_length = 0
        clustered = {}
        for i, hit in enumerate(self.fw_hits):
            pair = (hit[0], hit[1])
            if clustered.get(pair) is None:
                cluster = []
                cluster_len = 0
                cluster.append(hit)
                cluster_len += self.k
                last_x = hit[0]
                last_y = hit[1]
                clustered[pair] = True

                for j in range(i + 1, len(self.fw_hits)):

                    next_hit = self.fw_hits[j]

                    if next_hit[0] > last_x + self.L:
                        break

                    if next_hit[1] == last_y or next_hit[1] == last_y + 1:
                        pair = (next_hit[0], next_hit[1])
                        dist = pair[0] - last_x
                        last_x = next_hit[0]
                        last_y = next_hit[1]
                        cluster.append(next_hit)
                        clustered[pair] = True
                        if dist < self.k:
                            cluster_len += self.k - dist
                        else:
                            cluster_len += self.k

                if cluster_len > best_length:
                    best_length = cluster_len
                    best_cluster = cluster

        return best_cluster, best_length

    def rc_diag_group(self):
        best_cluster = None
        best_length = 0
        clustered = {}
        for i, hit in enumerate(self.rc_hits):
            pair = (hit[0], hit[1])
            if clustered.get(pair) is None:
                cluster = []
                cluster_len = 0
                cluster.append(hit)
                cluster_len += self.k
                last_x = hit[0]
                last_y = hit[1]
                clustered[pair] = True

                for j in range(i + 1, len(self.rc_hits)):
                    next_hit = self.rc_hits[j]
                    if next_hit[0] > last_x + self.L:
                        break

                    if next_hit[1] == last_y or next_hit[1] == last_y - 1:
                        pair = (next_hit[0], next_hit[1])
                        dist = pair[0] - last_x
                        last_x = next_hit[0]
                        last_y = next_hit[1]
                        cluster.append(next_hit)
                        clustered[pair] = True
                        if dist < self.k:
                            cluster_len += self.k - dist
                        else:
                            cluster_len += self.k

                if cluster_len > best_length:
                    best_length = cluster_len
                    best_cluster = cluster

        return best_cluster, best_length

    def cluster_hits(self, size_threshold = 5, group_hit = 1.0):
        self.chain_align = []
        query_len = self.length
        target_len = self.target_length

        align, length = self.fw_diag_group()
        self.fw_chain = align

        if align:
            # find left side extension length
            if align[0][0] < align[0][1] * self.L:
                left_extend = align[0][0]
            else:
                left_extend = align[0][1] * self.L

            if query_len - align[-1][0] < target_len - align[-1][1] * self.L:
                right_extend = query_len - align[-1][0]
            else:
                right_extend = target_len - align[-1][1] * self.L

            # middle span
            middle_extend = (abs(align[-1][0] - align[0][0]))

            extend = left_extend + middle_extend + right_extend

            if extend / float(group_hit * self.L) <= float(length) / self.k and length > size_threshold * self.k:
                if length > len(self.chain_align):
                    self.chain_align = align
                    self.is_forward = True

        align, length = self.rc_diag_group()
        self.rc_chain = align

        if align:
            # find left side extension length
            if align[0][0] < target_len - align[0][1] * self.k:
                left_extend = align[0][0]
            else:
                left_extend = target_len - align[0][1] * self.k

            if query_len - align[-1][0] < align[-1][1] * self.k:
                right_extend = query_len - align[-1][0]
            else:
                right_extend = align[-1][1] * self.k

            # middle span
            middle_extend = abs(align[-1][0] - align[0][0])

            extend = left_extend + middle_extend + right_extend

            if extend / float(group_hit * self.L) <= float(
                    length) / self.k and length > size_threshold * self.k:
                if length > len(self.chain_align):
                    self.chain_align = align
                    self.is_forward = False

        if len(self.chain_align) != 0:
            self.aligned = True
        else:
            self.aligned = False

    def plot(self):
        plt.figure()

        for hit in self.fw_hits:
            # hit[0] x coordinate, hit[3] list of y coordinate
            x = [hit[0]] * len(hit[3])
            plt.scatter(x, hit[3])


        for aligned_hit in self.fw_chain:
            x = [aligned_hit[0]] * len(aligned_hit[3])
            plt.scatter(x, aligned_hit[3], edgecolors="black", linewidths=2)


        plt.figure()

        for hit in self.rc_hits:
            # hit[0] x coordinate, hit[3] list of y coordinate
            x = [hit[0]] * len(hit[3])
            plt.scatter(x, hit[3])



        for aligned_hit in self.rc_chain:
            x = [aligned_hit[0]] * len(aligned_hit[3])
            plt.scatter(x, aligned_hit[3], edgecolors="black", linewidths=2)

        plt.show()


if __name__ == '__main__':
    # record1 = SeqIO.read("D:/Data/20170429/large_9mer_5_FP/FP_pair4_query.fasta", "fasta")
    # record2 = SeqIO.read("D:/Data/20170429/large_9mer_5_FP/FP_pair4_target.fasta", "fasta")
    record1 = SeqIO.read("D:/Data/20170614/FP_query_100.fasta", "fasta")
    record2 = SeqIO.read("D:/Data/20170614/FP_target_100.fasta", "fasta")
    test_filter = PseudoBloomFilter.PseudoBloomFilter(record2, 9, 135)
    print test_filter.L
    test_filter.generate_filter()
    test_query = QuerySeq(record1)
    test_query.check_kmer(test_filter)
    # print(test_query.fw_hits)
    # print(test_query.rc_hits)
    test_query.cluster_hits()
    print test_query.chain_align
    print test_query.aligned
    test_query.plot()

