import pickle
from Bio import SeqIO

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

blasr_large = {}
blasr_medium = {}
blasr_small = {}

# build seq dict so that we can based on index in output of daligener to find seq_id
query_seq = {}
query_records = SeqIO.parse("D:/Data/20170312/sensitivety_query.fasta","fasta")
for i,record in enumerate(query_records):
    query_seq[i+1] = record.id

large_seq = {}
large_records = SeqIO.parse("D:/Data/20170312/sensitivety_large_overlap.fasta","fasta")
for i,record in enumerate(large_records):
    large_seq[i+1] = record.id

medium_seq = {}
medium_records = SeqIO.parse("D:/Data/20170312/sensitivety_medium_overlap.fasta","fasta")
for i,record in enumerate(medium_records):
    medium_seq[i+1] = record.id

small_seq = {}
small_records = SeqIO.parse("D:/Data/20170312/sensitivety_small_overlap.fasta","fasta")
for i,record in enumerate(small_records):
    small_seq[i+1] = record.id

overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query.large_overlap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="q":
            #prin
            query_index = int(line[0:4])
            target_index = int(line[4:10])
            query_id = query_seq[query_index]
            target_id = large_seq[target_index]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id), False) is False:
                    #check is this pair tested before
                    num_found += 1
                    found[(query_id, target_id)] = True

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
print "sensitivity", float(true_align) / 610
print "FPR", (num_found - true_align) / float(42497- 610)

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query.medium_overlap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="q":
            #prin
            query_index = int(line[0:4])
            target_index = int(line[4:10])
            query_id = query_seq[query_index]
            target_id = medium_seq[target_index]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id), False) is False:
                    #check is this pair tested before
                    num_found += 1
                    found[(query_id, target_id)] = True

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
print "sensitivity", float(true_align) / 588
print "FPR", (num_found - true_align) / float(42499- 588)

num_found = 0
true_align = 0
found = {}
with open("D:/Data/20170406/query.small_overlap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="q":
            #prin
            query_index = int(line[0:4])
            target_index = int(line[4:10])
            query_id = query_seq[query_index]
            target_id = small_seq[target_index]
            if query_id!=target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id), False) is False:
                    #check is this pair tested before
                    num_found += 1
                    found[(query_id, target_id)] = True

                    if target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                            overlap_dict[query_id][2]:

                        true_align += 1

print num_found
print "sensitivity", float(true_align) / 558
print "FPR", (num_found - true_align) / float(42498- 558)