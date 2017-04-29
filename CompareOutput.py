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

# we start to build list of minimap output
minimap_output = []
found = {} # make sure no duplicate
with open("D:/Data/20170406/query_small_minimap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[5]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id),
                                                                                  False) is False:
                    minimap_output.append((query_id, target_id))
                    found[(query_id, target_id)] = True

# we build list of group hit output
group_output = []
with open("D:/Data/20170429/small_9mer_0.85_0.12_5_align_found.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            if query_id!=target_id:
                group_output.append((target_id, query_id))

print len(minimap_output)
print len(group_output)

minimap_extra = list(set(minimap_output) - set(group_output))
group_extra = list(set(group_output) - set(minimap_output))

TP_group_better = []
FP_group_worse = []
for pair in group_extra:
    query_id = pair[0]
    target_id = pair[1]
    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
            overlap_dict[query_id][2]:
        TP_group_better.append(pair)
    else:
        FP_group_worse.append(pair)


TP_minimap_better = []
FP_minimap_worse = []
for pair in minimap_extra:
    query_id = pair[0]
    target_id = pair[1]
    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
            overlap_dict[query_id][2]:
        TP_minimap_better.append(pair)
    else:
        FP_minimap_worse.append(pair)

with open("D:/Data/20170429/small_compare_with_minimap_9mer_0.85_0.12_5.out", "w") as fout:
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
