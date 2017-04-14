from Bio import SeqIO

overlap_pair = {}
expect = 0
with open("D:/Data/20170409/Ecoli_15X_to_ref_intersect.bed") as f_in:
    for line in f_in:
        sp = line.strip().split()
        query_id = sp[3]
        target_id = sp[9]
        if query_id != target_id:
            if overlap_pair.get((query_id, target_id), False) is False and overlap_pair.get((target_id, query_id), False) is False:
                overlap_pair[(query_id, target_id)] = True
                expect += 1

print expect

size = 0
for rec in SeqIO.parse("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq"):
    size += 1

total_pairs = size*(size - 1)/2

num_found = 0
true_align = 0
found = {}
true_align_pair = {}
false_align_pair = {}
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

                    if overlap_pair.get((query_id, target_id), False) is True or overlap_pair.get((target_id, query_id), False) is True:

                        # print target_id, query_id
                        true_align += 1
                        true_align_pair[(query_id, target_id)] = True
                    else:
                        false_align_pair[(query_id, target_id)] = True

print num_found
print true_align
print "sensitivity", float(true_align) / expect
print "FPR", (num_found - true_align) / float(total_pairs- expect)

with open("D:/Data/20170414/minimap_15X_k13w5_missing_align.out", "w") as fout:
    for pair in overlap_pair.keys():
        query = pair[0]
        target = pair[1]
        if true_align_pair.get((query, target), False) is False and true_align_pair.get((target, query),False) is False:
            print >> fout, pair[0] + "\t" + pair[1]

with open("D:/Data/20170414/minimap_15X_k13w5_fp_align.out", "w") as fout:
    for pair in false_align_pair:
        query = pair[0]
        target = pair[1]
        print >> fout, pair[0] + "\t" + pair[1]