from Bio import SeqIO
import pickle
from sklearn import metrics
import matplotlib.pyplot as plt


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

num_found = 0
true_align = 0
found = {}

with open("D:/Data/20170807/yass_9mer_small_overlap.out") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0] != "#" and line[0] != "(" and len(line) > 8:
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            bit_score = float(line_sp[11])
            if query_id != target_id:
                if found.get((query_id, target_id), False) is False and found.get((target_id, query_id),
                                                                                  False) is False:
                    # check is this pair tested before
                    num_found += 1
                    found[(query_id, target_id)] = bit_score

                    target_dict_list = overlap_dict.get(query_id)
                    if target_dict_list != None:

                        if target_id in target_dict_list[0] or target_id in target_dict_list[1] or target_id in \
                                target_dict_list[2]:
                            true_align += 1
                if found.get((query_id, target_id), False):
                    old_score = found.get((query_id, target_id), False)
                    if bit_score > old_score:
                        found[(query_id, target_id)] = bit_score

                if found.get((target_id, query_id), False):
                    old_score = found.get((target_id, query_id), False)
                    if bit_score > old_score:
                        found[(target_id, query_id)] = bit_score

print num_found
expect = 558
total_pairs = 42498
sensitivity = float(true_align) / expect
accuracy = float(true_align) / num_found
print total_pairs
print "sensitivity", sensitivity
print "FPR", (num_found - true_align) / float(total_pairs - expect)
print "accuracy", accuracy
print "F1", 2 * (accuracy * sensitivity) / (accuracy + sensitivity)

label = []
score = []
score_label = []
for pair in found.keys():
    if pair[0] != pair[1]:
        target_dict_list = overlap_dict.get(pair[0])
        if target_dict_list != None:

            if pair[1] in target_dict_list[0] or pair[1] in target_dict_list[1] or pair[1] in \
                    target_dict_list[2]:
                label.append(1)
            else:
                label.append(0)
        else:
            label.append(0)

        score.append(found[pair])
        score_label.append((score[-1], label[-1]))

fpr, tpr, thresholds = metrics.roc_curve(label, score, pos_label=1)
plt.plot(fpr, tpr)
plt.show()
score_label.sort(reverse=True)
print score_label
sensitvity_list = []
FPR_list = []
accuracy_list = []
previous_threshold = 10000.0
for i in range(len(score_label)):
    current_threshold = score_label[i][0]
    if current_threshold < previous_threshold:
        previous_threshold = current_threshold

        for j in range(len(score_label)):
            if score_label[j][0] < current_threshold:
                break

        true_align = sum([pair[1] for pair in score_label[0:j]])
        num_found = j
        sensitivity = float(true_align) / expect
        FPR = (num_found - true_align) / float(total_pairs - expect)
        accuracy_minus = 1 - float(true_align) / num_found
        sensitvity_list.append(sensitivity)
        FPR_list.append(FPR)
        accuracy_list.append(accuracy_minus)

plt.plot(FPR_list, sensitvity_list)
# plt.scatter(accuracy_list, sensitvity_list)
"""
plt.scatter(0.000861624291465, 0.953216374269)
plt.scatter(0.000873151037839,0.953216374269)
plt.scatter(0.00113826620444,0.988304093567)
plt.scatter(0.00107198741279,0.932748538012, label = "minimap")
"""
sen_yasschain = [
    0.433691756272,
    0.587813620072,
    0.743727598566,
    0.806451612903,
    0.836917562724,
    0.858422939068,
    0.874551971326,
    0.89247311828,
    0.905017921147,
    0.913978494624,
    0.921146953405
]
FPR_yasschain = [
    2.38435860753e-05,
    4.76871721507e-05,
    0.000143061516452,
    0.000166905102527,
    0.000190748688603,
    0.000238435860753,
    0.000309966618979,
    0.00035765379113,
    0.000834525512637,
    0.0015498330949,
    0.00448259418216
]
plt.scatter(FPR_yasschain,sen_yasschain)
sen_minimap = [
0.813620071685,
0.632616487455,
0.672043010753,
0.71146953405,
0.754480286738,
0.836917562724,
0.860215053763
]
FPR_minimap = [
0.000596089651884,
0.000214592274678,
0.000381497377206,
0.00035765379113,
0.00035765379113,
0.00100143061516,
0.0140677157845

]
plt.scatter(FPR_minimap, sen_minimap, label="minimap")
plt.legend()
plt.show()
