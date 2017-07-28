from bintrees import *
import ProbFunc


class GroupHit(object):
    def __init__(self, group_line):
        sp = group_line.strip().split("\t")
        self.query = sp[0]
        self.target = sp[1]
        self.direction = sp[2]

        # L and I is for DP chaining
        self.I = []
        # self.L = FastRBTree()


        groups = []
        for group in sp[4:]:
            hits = []
            group_sp = group.strip(",").split(",")
            for hit in group_sp:
                hit_sp = hit.split(" ")

                hits.append((int(hit_sp[0]) , int(hit_sp[1]) + int(hit_sp[0]), int(hit_sp[2])))

            self.I.append((hits[0][0] - hits[0][2], len(groups), 0))
            self.I.append((hits[-1][0], len(groups), -1))
            # self.L.insert((hits[-1][1], len(groups)))
            groups.append(hits)

        self.groups = groups

    def chain_groups(self ):

        r = len(self.I)
        L = FastRBTree()
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
                try:
                    j_key = L.floor_key(l_k - 1)

                    v_j = sum([x[2] for x in self.groups[k]])
                    j = L[j_key][1]

                    prev_score = V[j]

                except KeyError:
                    v_j = sum([x[2] for x in self.groups[k]])
                    j = -1
                    prev_score = 0

                V[k] = prev_score + v_j
                back_track[k] = j

            # is a end point
            else:
                #print L
                k = self.I[i][1]
                h_k = self.groups[k][-1][1]

                try:
                    j_key = L.ceiling_key(h_k)
                    j = L[j_key][1]
                    V_j = L[j_key][0]


                    if V[k] > V_j:
                        L.insert(h_k, (V[k], k))

                    #L.insert(h_k, (V[k], k))
                except KeyError:
                    # print len(L)
                    if len(L) == 0 or L.max_item()[1][0] < V[k]:
                        L[h_k] = (V[k], k)
                    # L.insert(h_k, (V[k], k))
                # max_item = L_tree.max_item()
                try:
                    j1_key = L.ceiling_key(h_k)
                    j1_score = L[j1_key]

                    # maybe mofify here'
                    #print "Vk", V[k]
                    #print L
                    while True:
                        prev_j1_key = j1_key
                        prev_j1_score = j1_score
                        try:
                            j1_key = L.succ_key(j1_key)
                            j1_score = L[j1_key]
                            if V[k] > j1_score[0]:
                                L.remove(j1_key)
                        except KeyError:
                            # print prev_j1_item
                            if V[k] > prev_j1_score[0]:
                                L.remove(prev_j1_key)
                            break

                    #print L
                except KeyError:
                    continue
        #print "DP finished"
        try:
            #print L
            max_item = L.max_item()
            score = max_item[1][0]

            current_j = max_item[1][1]
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


if __name__ == '__main__':
    file = "D:/Data/20170727/test.out"
    with open(file) as f:
        lines = f.readlines()

    for line in lines:
        group_hit = GroupHit(line)
        # print group_hit.groups
        print group_hit.chain_groups()
