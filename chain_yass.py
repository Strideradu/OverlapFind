from sortedcontainers import SortedDict
import ProbFunc
from matplotlib import pyplot as plt


def find_floor_key(d, key):
    try:
        i_loc = d.index(key)
        if i_loc == 0:
            return None, None, None
        else:
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]

    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == 0:
            del d[key]
            return None, None, None
        else:
            del d[key]
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]


def find_not_smaller_key(d, key):
    try:
        i_loc = d.index(key)
        return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]
    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == len(d) - 1:
            del d[key]
            return None, None, None
        else:
            # print i_loc
            # print len(d)
            del d[key]
            return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]


def find_not_bigger_key(d, key):
    try:
        i_loc = d.index(key)
        return i_loc, d.iloc[i_loc], d[d.iloc[i_loc]]
    except ValueError:
        # just insert something
        d[key] = 0
        i_loc = d.index(key)
        if i_loc == 0:
            del d[key]
            return None, None, None
        else:
            # print i_loc
            # print len(d)
            del d[key]
            return i_loc - 1, d.iloc[i_loc - 1], d[d.iloc[i_loc - 1]]


class GroupHit(object):
    def __init__(self, group_line):
        sp = group_line.strip().split("\t")
        self.target = sp[0]
        self.target_len = int(sp[1])
        self.query = sp[2]
        self.query_len = int(sp[3])
        self.aligned = False
        self.chain_align = None
        if sp[4] == 0:
            self.forward = True

        else:
            self.forward = False

        # L and I is for DP chaining
        self.I = []

        groups = []
        for group in sp[5:]:
            hits = []
            group_sp = group.strip(",").split(",")
            # print group_sp
            if len(group_sp) > 1 or int(group_sp[0].split()[2]) > 15:
                for hit in group_sp:
                    hit_sp = hit.split(" ")

                    if self.forward:
                        hits.append((int(hit_sp[0]), int(hit_sp[1]) + int(hit_sp[0]), int(hit_sp[2])))

                    else:
                        hits.append((int(hit_sp[0]), int(hit_sp[1]) - int(hit_sp[0]), int(hit_sp[2])))

                self.I.append((hits[0][0] - hits[0][2], len(groups), 0))
                self.I.append((hits[-1][0], len(groups), -1))
                # self.L.insert((hits[-1][1], len(groups)))
                groups.append(hits)

        self.groups = groups

    def fw_chain_groups(self):

        r = len(self.I)
        # L = FastRBTree()
        L_dict = SortedDict()
        V = [0] * len(self.groups)
        back_track = [-1] * len(self.groups)

        self.I.sort()
        # print self.I
        ## print self.groups

        for i in range(r):
            # I_list[i] is a start point, noticed we go through the chain from botton to top
            if self.I[i][2] == 0:
                k = self.I[i][1]

                l_k = self.groups[k][0][1] - self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0] - self.groups[k][0][2]
                end_x = self.groups[k][-1][0]
                end_y = self.groups[k][-1][1]
                diagonal = start_y - start_x
                # find largest h_j strictly smaller than l_k and also not off diagonal

                j_index, j_key, j_value = find_floor_key(L_dict, l_k)
                # j_key = L.floor_key(l_k - 1)
                if j_index != None:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = j_value[1]

                    prev_score = V[j]

                else:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = -1
                    prev_score = 0

                V[k] = prev_score + v_j
                back_track[k] = j

            # is a end point
            else:
                # print L
                k = self.I[i][1]
                h_k = self.groups[k][-1][1]

                j_index, j_key, j_value = find_not_smaller_key(L_dict, h_k)
                # j_key = L.ceiling_key(h_k)
                if j_index != None:
                    j = j_value[1]
                    V_j = j_value[0]

                    if V[k] > V_j:
                        L_dict[h_k] = (V[k], k)
                        # L.insert(h_k, (V[k], k))

                # L.insert(h_k, (V[k], k))
                else:
                    # print len(L)
                    if len(L_dict) == 0 or L_dict.peekitem()[1][0] < V[k]:
                        L_dict[h_k] = (V[k], k)
                """
                if len(L) == 0 or L.max_item()[1][0] < V[k]:
                    L[h_k] = (V[k], k)
                """
                # L.insert(h_k, (V[k], k))

                j1_index, j1_key, j1_value = find_not_smaller_key(L_dict, h_k)
                # maybe mofify here'
                # print "Vk", V[k]
                # print L
                if j1_index != None:
                    j1_index += 1
                    while j1_index < len(L_dict):
                        prev_j1_key = j1_key
                        prev_j1_value = j1_value
                        try:

                            j1_key = L_dict.iloc[j1_index]
                            j1_value = L_dict[j1_key]
                            if V[k] > j1_value[0]:
                                del L_dict[j1_key]
                            else:
                                j1_index += 1
                        except KeyError:
                            # print prev_j1_item
                            if V[k] > prev_j1_value[0]:
                                del L_dict[prev_j1_key]
                            break

        # print "DP finished"
        try:
            # print L
            max_item = L_dict.peekitem(index=-1)
            max_value = max_item[1]
            score = max_value[0]

            current_j = max_value[1]
            # backtrack
            chain_index = []
            chain_index.append(current_j)

            while True:
                prev_j = back_track[current_j]
                if prev_j == -1:
                    break
                else:
                    current_j = prev_j
                    chain_index.append(current_j)

            optimal_chain = []
            length = 0
            for i in chain_index[::-1]:
                optimal_chain.extend(self.groups[i])
                length += sum([x[2] for x in self.groups[i]])

        except ValueError:
            optimal_chain = None
            length = 0

        return optimal_chain, length

    def rc_chain_groups(self):

        r = len(self.I)
        # L = FastRBTree()
        L_dict = SortedDict()
        V = [0] * len(self.groups)
        back_track = [-1] * len(self.groups)

        self.I.sort()
        # print self.I
        ## print self.groups

        for i in range(r):
            # I_list[i] is a start point, noticed we go through the chain from botton to top
            if self.I[i][2] == 0:
                k = self.I[i][1]

                l_k = self.groups[k][0][1] - self.groups[k][0][2]
                start_y = l_k
                start_x = self.I[i][0] - self.groups[k][0][2]
                end_x = self.groups[k][-1][0]
                end_y = self.groups[k][-1][1]
                diagonal = start_y + start_x
                # find largest h_j strictly smaller than l_k and also not off diagonal

                j_index, j_key, j_value = find_not_smaller_key(L_dict, l_k + 1)
                # j_key = L.floor_key(l_k - 1)
                if j_index != None:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = j_value[1]

                    prev_score = V[j]

                else:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = -1
                    prev_score = 0

                V[k] = prev_score + v_j
                back_track[k] = j

            # is a end point
            else:
                # print L
                k = self.I[i][1]
                h_k = self.groups[k][-1][1]

                j_index, j_key, j_value = find_not_bigger_key(L_dict, h_k)
                # j_key = L.ceiling_key(h_k)
                if j_index != None:
                    j = j_value[1]
                    V_j = j_value[0]

                    if V[k] > V_j:
                        L_dict[h_k] = (V[k], k)
                        # L.insert(h_k, (V[k], k))

                # L.insert(h_k, (V[k], k))
                else:
                    # print len(L)
                    if len(L_dict) == 0 or L_dict.peekitem(0)[1][0] < V[k]:
                        L_dict[h_k] = (V[k], k)
                """
                if len(L) == 0 or L.max_item()[1][0] < V[k]:
                    L[h_k] = (V[k], k)
                """
                # L.insert(h_k, (V[k], k))

                j1_index, j1_key, j1_value = find_not_bigger_key(L_dict, h_k)
                # maybe mofify here'
                # print "Vk", V[k]
                # print L
                if j1_index != None:
                    j1_index -= 1
                    while j1_index >= 0:
                        prev_j1_key = j1_key
                        prev_j1_value = j1_value
                        try:
                            # print j1_index
                            j1_key = L_dict.iloc[j1_index]
                            j1_value = L_dict[j1_key]
                            if V[k] > j1_value[0]:
                                del L_dict[j1_key]

                            j1_index -= 1
                        except KeyError:
                            # print prev_j1_item
                            if V[k] > prev_j1_value[0]:
                                del L_dict[prev_j1_key]
                            break

        # print "DP finished"
        try:
            # print L
            max_item = L_dict.peekitem(index=0)
            max_value = max_item[1]
            score = max_value[0]

            current_j = max_value[1]
            # backtrack
            chain_index = []
            chain_index.append(current_j)

            while True:
                prev_j = back_track[current_j]
                if prev_j == -1:
                    break
                else:
                    current_j = prev_j
                    chain_index.append(current_j)

            optimal_chain = []
            length = 0
            for i in chain_index[::-1]:
                optimal_chain.extend(self.groups[i])
                length += sum([x[2] for x in self.groups[i]])

        except ValueError:
            optimal_chain = None
            length = 0

        return optimal_chain, length

    def chain_groups(self, accuracy=0.8, group_distance=None, rechain_threshold=5, span_coefficient=1.0):
        if group_distance == None:
            group_distance = ProbFunc.statistical_bound_of_waiting_time(accuracy, 9)

        if len(self.groups) != 0:
            extend = None
            if self.forward:
                align, length = self.fw_chain_groups()

                if align:
                    # find left side extension length
                    if align[0][0] < align[0][1]:
                        left_extend = align[0][0] - align[0][2]
                    else:
                        left_extend = align[0][1] - align[0][2]

                    if self.query_len - align[-1][0] < self.target_len - align[-1][1]:
                        right_extend = self.query_len - align[-1][0]
                    else:
                        right_extend = self.target_len - align[-1][1]

                    # middle span
                    middle_extend = 0.5 * (abs(align[-1][0] - align[0][0] + align[0][2])) + 0.5 * abs(
                        align[-1][1] - align[0][1] + align[0][2])

                    extend = left_extend + middle_extend + right_extend

                """
                extend = 2 * group_distance + 0.5 * (
                    align[-1][0] - align[0][0] + align[0][2] + align[-1][1] - align[0][1] + align[0][2])
                """

            else:
                align, length = self.rc_chain_groups()

                if align:
                    # find left side extension length
                    if align[0][0] < self.target_len - align[0][1]:
                        left_extend = align[0][0] - align[0][2]
                    else:
                        left_extend = self.target_len - align[0][1]

                    if self.query_len - align[-1][0] < align[-1][1]:
                        right_extend = self.query_len - align[-1][0]
                    else:
                        right_extend = align[-1][1]

                    # middle span
                    middle_extend = 0.5 * (abs(align[-1][0] - align[0][0] + align[0][2])) + 0.5 * abs(
                        align[0][1] - align[-1][1] + align[0][2])

                    extend = left_extend + middle_extend + right_extend
                    #print left_extend
                    #print right_extend
                    #print self.query_len
                    #print self.target_len

                """
                extend = 2 * group_distance + 0.5 * (
                    align[-1][0] - align[0][0] + align[0][2] + align[0][1] - align[-1][1] + align[0][2])
                """
            if extend:
                #print align
                #print extend
                if extend / float(span_coefficient * group_distance) <= float(
                        length) / 9 and length > rechain_threshold * 9:
                    self.chain_align = align
                    self.aligned = True

    def plot_hits(self):
        plt.figure()

        if self.aligned:
            coordinates = map(list, zip(*self.chain_align))
            plt.scatter(coordinates[0], coordinates[1], s = 40,edgecolors="black", linewidths=5)
        for group in self.groups:
            # print chain
            chain_coor = map(list, zip(*group))
            plt.scatter(chain_coor[0], chain_coor[1], s=40)



        plt.show()


if __name__ == '__main__':
    file = "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170731/query_all_0p85_0p12.out"
    # file = "D:/Data/20170727/query_all_0p85_0p12_first10.out"
    #file = "C:/Users/Nan/Documents/GitHub/yass/cmake-build-debug/fp_002.out"
    L = ProbFunc.statistical_bound_of_waiting_time(0.85, 9)

    with open(file) as f:
        lines = f.readlines()

    for line in lines:
        group_hit = GroupHit(line)
        # print group_hit.groups

        group_hit.chain_groups(accuracy=0.85, group_distance=L, rechain_threshold=3, span_coefficient=1.0)
        print group_hit.chain_align
        if group_hit.aligned:
            output_str = group_hit.query + "\t" + group_hit.target + "\t" + str(group_hit.aligned)
            print(output_str)

        group_hit.plot_hits()
