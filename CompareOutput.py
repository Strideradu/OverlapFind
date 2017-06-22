import pickle

# Trying to find the difference of output from minimap and our group hit strategies
# we need 4 sets: Group Hit found true align, but minimap fails
#                 minimap found true align, but group hit fails
#                 minimap has FP but group hit not
#                 group hit has FP but minimap not

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")

group_output = []
found = {}
with open("D:/Data/20170622/query_all_target_9_0.85_3_new.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="#" and line[0]!="(":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id),
                                                                                  False) is False:

                    group_output.append((query_id, target_id))
                    found[(query_id, target_id)] = True

# we start to build list of minimap output
minimap_output = []
minimap_found = {} # make sure no duplicate
with open("D:/Data/20170613/query_all_minimap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[5]
            if query_id!=target_id:
                if minimap_found.get((query_id, target_id), False) is False and minimap_found.get((target_id, query_id),
                                                                                  False) is False:
                    if found.get((target_id, query_id), False) is False:
                        minimap_output.append((query_id, target_id))
                    else:
                        minimap_output.append((target_id, query_id))
                    minimap_found[(query_id, target_id)] = True

# we build list of group hit output


print len(minimap_output)
print len(group_output)

minimap_extra = list(set(minimap_output) - set(group_output))
group_extra = list(set(group_output) - set(minimap_output))

TP_group_better = []
FP_group_worse = []
for pair in group_extra:
    query_id = pair[0]
    target_id = pair[1]
    if overlap_dict.get(query_id) and (target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
            overlap_dict[query_id][2]):
        TP_group_better.append(pair)
    else:
        FP_group_worse.append(pair)


TP_minimap_better = []
FP_minimap_worse = []
for pair in minimap_extra:
    query_id = pair[0]
    target_id = pair[1]
    if overlap_dict.get(query_id) and (target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
            overlap_dict[query_id][2]):
        TP_minimap_better.append(pair)
    else:
        FP_minimap_worse.append(pair)

with open("D:/Data/20170622/compare_with_minimap_9mer_0.85_3_corrected_extend.out", "w") as fout:
    print >> fout, "# True Align only found by Group Hit:"
    for pair in TP_group_better:
        print >> fout, pair[0] + "\t" + pair[1]

    print >> fout, "# True Align only found by Minimap:"
    for pair in TP_minimap_better:
        print >> fout, pair[0] + "\t" + pair[1]

    print >> fout, "# False Positive only found by Group Hit:"
    for pair in FP_group_worse:
        print >> fout, pair[0] + "\t" + pair[1]

    print >> fout, "# False Positive only found by Minimap:"
    for pair in FP_minimap_worse:
        print >> fout, pair[0] + "\t" + pair[1]
