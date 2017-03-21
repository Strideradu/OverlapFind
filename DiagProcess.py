"""
A simple class to find all kmers and process to generate diagonal plot
didn't deal the exactly same case
"""
from matplotlib import pyplot as plt
from Bio import SeqIO
from QualitySeq import QualitySeq
import ProbFunc

def reverse_com(string):
    rev_com = {"A":"T", "C":"G", "T":"A", "G":"C"}
    result = ""
    for i in string:
        result = rev_com[i] + result

    return result


class DiagProcess(object):
    def __init__(self, seq1, seq2):
        self.query = seq1
        self.target = seq2
        self.fw_points = None
        self.rc_points = None
        self.k = 0
        self.fw_chain = None
        self.rc_chain = None
        self.chain_align = []
        self.is_forward = None
        self.fw_seeds = 0
        self.rc_seeds = 0
        self.aligned = None


    def diag_points(self, k):
        self.k = k

        dict1 = self.query.generate_kmer_pos(k)
        dict2 = self.target.generate_kmer_pos(k)

        self.fw_points = []
        self.rc_points = []

        allkey = set().union(dict1, dict2)

        for key in allkey:
            # forward case
            if not(dict2.get(key) is None) and not(dict1.get(key) is None):
                query_list = dict1[key]
                target_list = dict2[key]

                for x_i in query_list:
                    for y_j in target_list:
                        point = (x_i[0], y_j[0], x_i[1] + y_j[1])
                        self.fw_points.append(point)

            # reverse complement case
            rc_key = reverse_com(key)
            if not(dict2.get(rc_key) is None) and not(dict1.get(key) is None):
                query_list = dict1[key]
                target_list = dict2[rc_key]

                for x_i in query_list:
                    for y_j in target_list:
                        point = (x_i[0] - k, y_j[0], x_i[1] + y_j[1])
                        self.rc_points.append(point)
            self.fw_points.sort()
            self.rc_points.sort()

    def diag_chain(self, accuracy = 0.75, gap = 0.2):
        """
        chaining the seeds
        :param accuracy: expected match rate
        :param gap: expected gap rate
        :return:
        """

        L = ProbFunc.statistical_bound_of_waiting_time(accuracy, self.k)
        delta = ProbFunc.statistical_bound_of_randomwalk(gap, L)

        # L = 64
        # delta = 6

        # process the forward read comparison
        # the list of chain we found, tuple of (seeds, l), seed is the list of all seed in the chain, l is the aligned lengthh.
        fw_chain = []
        chained = {}

        i = 0
        while i < len(self.fw_points):
            if chained.get(self.fw_points[i], False) is False:
                # print i
                chain = []
                l = self.k
                last_i = i
                last_x = self.fw_points[i][0]
                last_y = self.fw_points[i][1]
                last_diagonal = last_y - last_x
                chain.append(self.fw_points[i])
                chained[self.fw_points[i]] = True
                for j in range(i + 1, len(self.fw_points)):
                    x = self.fw_points[j][0]
                    y = self.fw_points[j][1]
                    diagonal =  y - x

                    delta_y = abs(y - last_y)
                    delta_x = abs(x - last_x)

                    if max(delta_y, delta_x)<=L and abs(diagonal - last_diagonal)<=delta:
                        if max(delta_y, delta_x)<= self.k:
                            l += max(delta_y, delta_x)
                        else:
                            l += self.k
                        last_x = x
                        last_y = y
                        last_diagonal = diagonal
                        chain.append(self.fw_points[j])
                        chained[self.fw_points[j]] = True


                    elif min(delta_y, delta_x)>L:
                        last_i = j
                        break

                if len(chain) > 1 and (abs(chain[-1][0] - chain[0][0]) > self.k or abs(chain[-1][1] - chain[0][1]) > self.k):
                    fw_chain.append((chain, l))
                    self.fw_seeds += len(chain)
                # print chain[0]

            i += 1

        self.fw_chain = fw_chain

        # process the reverse chaining

        rc_chain = []
        chained = {}

        i = 0
        while i < len(self.rc_points):
            if chained.get(self.rc_points[i], False) is False:
                # print i
                chain = []
                last_i = i
                l = self.k
                last_x = self.rc_points[i][0]
                last_y = self.rc_points[i][1]
                last_diagonal = last_y + last_x
                chain.append(self.rc_points[i])
                chained[self.rc_points[i]] = True
                for j in range(i + 1, len(self.rc_points)):
                    x = self.rc_points[j][0]
                    y = self.rc_points[j][1]
                    diagonal =  y + x

                    delta_y = abs(y - last_y)
                    delta_x = abs(x - last_x)

                    if max(delta_y, delta_x)<=L and abs(diagonal - last_diagonal)<=delta:
                        if max(delta_y, delta_x)<= self.k:
                            l += max(delta_y, delta_x)
                        else:
                            l += self.k
                        last_x = x
                        last_y = y
                        last_diagonal = diagonal
                        chain.append(self.rc_points[j])
                        chained[self.rc_points[j]] = True

                    elif min(delta_y, delta_x)>L:
                        # last_i = j
                        break

                if len(chain) > 1 and (abs(chain[-1][0] - chain[0][0]) > self.k and abs(chain[-1][1] - chain[0][1]) > self.k):
                    rc_chain.append((chain, l))
                    self.rc_seeds += len(chain)
                # print chain[0]

            i += 1

        self.rc_chain = rc_chain

    def rechain(self, gap = 0.2, chain_threshold=5):
        """
        connect the chain of the seeds from chaining step to generate final long chain
        :return:
        """
        # print "start"
        # connect the forward chains
        connected = {}

        i = 0
        # print len(self.fw_chain)
        seed_num = self.fw_seeds
        previous_seed_num = 0
        while i < len(self.fw_chain):
            #print i
            if connected.get(i, False) is False:
                align = []
                length  = 0
                last_end_x = self.fw_chain[i][0][-1][0]
                last_end_y = self.fw_chain[i][0][-1][1]

                last_end_diagonal = last_end_y - last_end_x
                align.extend(self.fw_chain[i][0])
                length += self.fw_chain[i][1]
                connected[i] = True
                for j in range(i + 1, len(self.fw_chain)):
                    #print "j", j
                    start_x = self.fw_chain[j][0][0][0]
                    start_y = self.fw_chain[j][0][0][1]
                    diagonal =  start_y - start_x
                    delta_y = abs(start_y - last_end_y)
                    delta_x = abs(start_x - last_end_x)
                    delta = max(delta_x, delta_y) * gap
                    # print delta

                    if abs(diagonal - last_end_diagonal) <= delta:
                        last_end_x = self.fw_chain[j][0][-1][0]
                        last_end_y = self.fw_chain[j][0][-1][1]
                        last_end_diagonal = last_end_y - last_end_x
                        align.extend(self.fw_chain[j][0])
                        length += self.fw_chain[j][1]
                        connected[j] = True

                if len(align) > len(self.chain_align):
                    # 9 mer with 0.75 threshold = 135

                    x_span = abs(align[-1][0] - align[0][0])
                    y_span = abs(align[-1][1] - align[0][1])
                    #print length / self.k
                    #print max(x_span, y_span) / 135

                    if max(x_span, y_span)/135 < 3 * length/self.k and length > chain_threshold*self.k:
                        self.chain_align = align
                        self.is_forward = True

            if seed_num - len(self.chain_align) - previous_seed_num < 0:
                break

            previous_seed_num += len(self.fw_chain[i][0])

            i += 1

        connected = {}

        i = 0
        seed_num = self.rc_seeds
        previous_seed_num = 0
        while i < len(self.rc_chain) - 1:
            if connected.get(i, False) is False:
                align = []
                length = 0
                last_end_x = self.rc_chain[i][0][-1][0]
                last_end_y = self.rc_chain[i][0][-1][1]

                last_end_diagonal = last_end_y + last_end_x
                align.extend(self.rc_chain[i][0])
                length += self.rc_chain[i][1]
                connected[i] = True
                for j in range(i + 1, len(self.rc_chain)):
                    start_x = self.rc_chain[j][0][0][0]
                    start_y = self.rc_chain[j][0][0][1]
                    diagonal =  start_y + start_x
                    delta_y = abs(start_y - last_end_y)
                    delta_x = abs(start_x - last_end_x)
                    delta = max(delta_x, delta_y) * gap

                    if abs(diagonal - last_end_diagonal) <= delta:
                        last_end_x = self.rc_chain[j][0][-1][0]
                        last_end_y = self.rc_chain[j][0][-1][1]
                        last_end_diagonal = last_end_y + last_end_x
                        align.extend(self.rc_chain[j][0])
                        length += self.rc_chain[j][1]
                        connected[j] = True

                if len(align) > len(self.chain_align):

                    x_span = abs(align[-1][0] - align[0][0])
                    y_span = abs(align[-1][1] - align[0][1])
                    #print length / self.k
                    #print max(x_span, y_span) / 135
                    # print x_span, y_span
                    if max(x_span, y_span) / 135 < 3 * length / self.k and length > chain_threshold*self.k:
                        self.chain_align = align
                        self.is_forward = False

            # if the remaning seed number is too small, that menas we already find longest align
            if seed_num - len(self.chain_align) - previous_seed_num < 0:
                break

            previous_seed_num += len(self.rc_chain[i][0])

            i += 1
        if len(self.chain_align) != 0:
            self.aligned = True
        else:
            self.aligned = False



    def diag_plot(self):
        plt.figure()
        coordinates = map(list, zip(*self.fw_points))

        plt.scatter(coordinates[0], coordinates[1], c=coordinates[2], cmap="Reds")
        plt.colorbar()
        for chain in self.fw_chain:
            #print chain
            chain_coor = map(list, zip(*chain[0]))
            plt.scatter(chain_coor[0], chain_coor[1], edgecolors="black", linewidths=2)
        plt.figure()
        coordinates = map(list, zip(*self.rc_points))
        plt.scatter(coordinates[0], coordinates[1],c = coordinates[2], cmap="Greens")
        plt.colorbar()
        for chain in self.rc_chain:
            #print chain
            chain_coor = map(list, zip(*chain[0]))
            plt.scatter(chain_coor[0], chain_coor[1], edgecolors="black", linewidths=2)
        plt.show()

if __name__ == '__main__':
    #record1 = SeqIO.read("D:/Data/20170213/unaligned_pair_3_1.fastq", "fastq")
    #record2 = SeqIO.read("D:/Data/20170213/unaligned_pair_3_2.fastq", "fastq")
    record1 = SeqIO.read("D:/Data/20170213/pair1_query.fastq", "fastq")
    record2 = SeqIO.read("D:/Data/20170213/pair1_target.fastq", "fastq")
    #record1 = SeqIO.read("D:/Data/20170321/Flase_Positive_Pair1_1.fastq", "fastq")
    #record2 = SeqIO.read("D:/Data/20170321/Flase_Positive_Pair1_4.fastq", "fastq")
    seq1 = QualitySeq(record1)
    seq2 = QualitySeq(record2)
    process = DiagProcess(seq1, seq2)
    process.diag_points(9)
    chians = process.diag_chain(0.75, 0.2)
    process.rechain(0.2)
    print process.chain_align
    print process.aligned
    process.diag_plot()
