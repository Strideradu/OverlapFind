"""
Query sequence process
"""
from Bio import SeqIO
import PseudoBloomFilter
from matplotlib import pyplot as plt
import math

gap_rate = 0.05

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

                    dist = next_hit[0] - cluster[0][0]
                    diag_off = min(dist * gap_rate, 2*self.L)
                    expect_y_max = math.ceil(float((cluster[0][1] -0.5) * self.L + dist + diag_off)/self.L)
                    expect_y_min = math.ceil(float((cluster[0][1] - 0.5) * self.L + dist - diag_off) / self.L)

                    if next_hit[1] >= expect_y_min and next_hit[1] <= expect_y_max:
                        pair = (next_hit[0], next_hit[1])

                        last_x = next_hit[0]
                        last_y = next_hit[1]
                        cluster.append(next_hit)
                        clustered[pair] = True
                        if dist < self.k:
                            cluster_len += dist
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

                    dist = next_hit[0] - cluster[0][0]
                    diag_off = min(dist * gap_rate, 2 * self.L)
                    expect_y_min = math.ceil(float((cluster[0][1] - 0.5) * self.L - dist - diag_off) / self.L)
                    expect_y_max = math.ceil(float((cluster[0][1] - 0.5) * self.L - dist + diag_off) / self.L)

                    if next_hit[1] >= expect_y_min and next_hit[1] <= expect_y_max:
                        pair = (next_hit[0], next_hit[1])
                        dist = pair[0] - last_x
                        last_x = next_hit[0]
                        last_y = next_hit[1]
                        cluster.append(next_hit)
                        clustered[pair] = True
                        if dist < self.k:
                            cluster_len += dist
                        else:
                            cluster_len += self.k

                if cluster_len > best_length:
                    best_length = cluster_len
                    best_cluster = cluster

        return best_cluster, best_length

    def cluster_hits(self, size_threshold = 5, debug= False, group_hit = 1.0):
        self.chain_align = []
        self.aligned = False
        query_len = self.length
        target_len = self.target_length

        align, length = self.fw_diag_group()
        if debug:
            print "fw_align_length:", length
            print align
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

            if debug:
                print "extend", extend

            x_extend = align[-1][0] - align[0][0]
            y_extend = (align[-1][1] - align[0][1] + 1) * self.L
            diag_off = int(x_extend * 0.2+1)

            if align[0][1] != align[-1][1]:
                if extend / float(2 * group_hit * self.L) <= float(length) / self.k and length > size_threshold * self.k:
                    if length > len(self.chain_align):
                        self.chain_align = align
                        self.is_forward = True

        align, length = self.rc_diag_group()
        if debug:
            print "rc_align_length:", length
            print align
        self.rc_chain = align

        if align:
            # find left side extension length
            if align[0][0] < target_len - align[0][1] * self.L:
                left_extend = align[0][0]
            else:
                left_extend = target_len - align[0][1] * self.L

            # print left_extend

            if query_len - align[-1][0] < align[-1][1] * self.L:
                right_extend = query_len - align[-1][0]
            else:
                right_extend = align[-1][1] * self.L

            # print right_extend

            # middle span
            middle_extend = abs(align[-1][0] - align[0][0])

            extend = left_extend + middle_extend + right_extend

            if debug:
                print "extend", extend

            x_extend = align[-1][0] - align[0][0]
            y_extend = (align[0][1] - align[-1][1]) * self.k
            diag_off = x_extend*0.2

            if align[0][1]!= align[-1][1]:
                if extend / float(2 * group_hit * self.L) <= float(
                        length) / self.k and length > size_threshold * self.k:
                    if length > len(self.chain_align):
                        self.chain_align = align
                        self.is_forward = False

        if len(self.chain_align) != 0:
            self.aligned = True


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
    # record1 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair1_query.fasta", "fasta")
    # record2 = SeqIO.read("D:/Data/20170429/large_9mer_5_missing/missing_pair1_target.fasta", "fasta")
    # record1 = SeqIO.read("D:/Data/20170622/9mer_FP/FP_query_002.fasta", "fasta")
    # record2 = SeqIO.read("D:/Data/20170622/9mer_FP/FP_target_002.fasta", "fasta")
    record1 = SeqIO.read("D:/Data/20170622/9mer_missing/missing_query_002.fasta", "fasta")
    record2 = SeqIO.read("D:/Data/20170622/9mer_missing/missing_target_002.fasta", "fasta")
    test_filter = PseudoBloomFilter.PseudoBloomFilter(record2, 9, 54)
    print test_filter.L
    test_filter.generate_filter()
    test_query = QuerySeq(record1)
    test_query.check_kmer(test_filter)
    # print(test_query.fw_hits)
    # print(test_query.rc_hits)
    test_query.cluster_hits(size_threshold = 3, debug = True)
    print test_query.chain_align
    print test_query.aligned
    test_query.plot()
