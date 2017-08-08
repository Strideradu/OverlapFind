import pickle
from Bio import SeqIO
# minimap output format, tab delimited file with the following column
# query name, length, 0-based start, end, strand, target name, length, start, end, the number of matching bases,
#  the number of co-linear minimizers in the match, the fraction of matching bases

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query_large_minimap.out") as f1:
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

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
sensitivity = float(true_align) / 610
accuracy = float(true_align)/num_found
print "sensitivity", sensitivity
print "FPR", (num_found - true_align) / float(42497- 610)
print "accuracy", accuracy
print "F1", 2*(accuracy*sensitivity)/(accuracy+sensitivity)

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query_medium_minimap.out") as f1:
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

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
sensitivity = float(true_align) / 588
accuracy = float(true_align)/num_found
print "sensitivity", sensitivity
print "FPR", (num_found - true_align) / float(42499- 588)
print "accuracy", accuracy
print "F1", 2*(accuracy*sensitivity)/(accuracy+sensitivity)

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query_small_minimap.out") as f1:
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

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
sensitivity = float(true_align) / 558
accuracy = float(true_align)/num_found
print "sensitivity", sensitivity
print "FPR", (num_found - true_align) / float(42498- 558)
print "accuracy", accuracy
print "F1", 2*(accuracy*sensitivity)/(accuracy+sensitivity)