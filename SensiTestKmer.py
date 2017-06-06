from Bio import SeqIO
import pickle


def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


size = 0
seq_id = []
for rec in SeqIO.parse("D:/Data/20170523/target_masked.fasta", "fasta"):
    seq_id.append(rec.id)
    size += 1

total_pairs = size * (size - 1) / 2

expect = 0
overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")
found_overlap = {}
for seq in overlap_dict.keys():
    if seq in seq_id:
        for target in (overlap_dict[seq][0] + overlap_dict[seq][1] + overlap_dict[seq][2]):
            if target in seq_id and found_overlap.get((seq, target)) is None and found_overlap.get(
                    (target, seq)) is None:
                expect += 1
                found_overlap[(seq, target)] = True
print expect

hit_thresholds = [0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50]

for hit_threshold in hit_thresholds:
    num_found = 0
    true_align = 0
    found = {}
    with open("D:/Data/20170603/compare_kmer_11mer.out") as f1:
        for line in f1:
            line = line.rstrip()
            if line != "" and line[0] != "#":
                line_sp = line.split("\t")
                query_id = line_sp[0]
                target_id = line_sp[1]
                num_hit = int(line_sp[2])
                if query_id != target_id:
                    if found.get((query_id, target_id), False) is False and found.get((target_id, query_id),
                                                                                      False) is False and num_hit > hit_threshold:
                        # check is this pair tested before
                        num_found += 1
                        found[(query_id, target_id)] = True

                        target_dict_list = overlap_dict.get(query_id)
                        if target_dict_list != None:

                            if target_id in target_dict_list[0] or target_id in target_dict_list[1] or target_id in \
                                    target_dict_list[2]:
                                true_align += 1
    print
    print "Hit Threshold:", hit_threshold
    print num_found
    sensitivity = float(true_align) / expect
    accuracy = float(true_align) / num_found
    print "sensitivity", sensitivity
    print "FPR", (num_found - true_align) / float(total_pairs - expect)
    print "accuracy", accuracy
    print "F1", 2 * (accuracy * sensitivity) / (accuracy + sensitivity)
    print "reduce", total_pairs/float(num_found)
