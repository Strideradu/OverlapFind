from Bio import SeqIO
import pickle

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

expect = 0
overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")
for seq in overlap_dict.keys():
    expect += len(overlap_dict[seq][0] + overlap_dict[seq][1] + overlap_dict[seq][2])

expect = expect/2.0
print expect

size = 0
for rec in SeqIO.parse("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq"):
    size += 1

total_pairs = size*(size - 1)/2

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170407/minimap_all_vs_all_15X_k13_w5.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[5]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id), False) is False:
                    #check is this pair tested before
                    num_found += 1
                    found[(query_id, target_id)] = True
                    try:

                        if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                                overlap_dict[query_id][2]:
                            # print target_id, query_id
                            true_align += 1

                    except KeyError:
                        continue

print num_found
print true_align
print "sensitivity", float(true_align) / expect
print "FPR", (num_found - true_align) / float(total_pairs- expect)