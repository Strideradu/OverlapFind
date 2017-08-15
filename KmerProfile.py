import matplotlib
matplotlib.use('Agg')
from Bio import SeqIO
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess
import pickle
from matplotlib import pyplot as plt
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("query", help="query fasta path", type=str)
parser.add_argument("target", help="target fasta path", type=str)
parser.add_argument("fig", help="fig path", type=str)
parser.add_argument("k", help="kmer", type=int)
parser.add_argument("max_k", help="max shared kmer in figure", type=int)
parser.add_argument("bin_size", help="number of bins", type=int)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

overlap_dict = load_obj("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/overlap.pkl")


# masked_fasta=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170324_dust_group_hit/filtered_15X_masked.fasta", "fasta")

query_fasta = SeqIO.parse(args.query, "fasta")

# target only extract the id, sequence should get from masked_fasta
target_fasta = SeqIO.parse(args.target, "fasta")
# SeqIO only generate an iterator so cannot iterate many times, generate a list first
target_list = list(target_fasta)
"""
for target_seq in target_fasta:
    target_list.append(target_seq.id)
"""
true = []
false = []
tested_pair = 0
for query_seq in query_fasta:
    seq_1 = QualitySeq(query_seq)
    query_id = query_seq.id
    # print query_seq.id

    for target_seq in target_list:
        target_id = target_seq.id
        if target_id != query_id:
            tested_pair += 1
            # target_masked = masked_fasta[target_id]
            seq_2 = QualitySeq(target_seq)
            process = DiagProcess(seq_1, seq_2)
            process.diag_points(args.k)
            if len(process.fw_points) > len(process.rc_points):
                shared_num = len(process.fw_points)
            else:
                shared_num = len(process.rc_points)

            # determine if they overlap each other

            if overlap_dict.get(query_id, False) and (target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][1] or target_id in \
                    overlap_dict[query_id][2]):
                true.append(shared_num)

            else:
                false.append(shared_num)

true.sort()
false.sort()
true_count = len(true)
false_count = len(false)
max_value = max(true[-1], false[-1])

sen_list = []
FPR_list = []

j1 = 0
j2 = 0
for i in range(max_value):
    while j1 < true_count:
        if true[j1] > i:
            j1 = j1 - 1
            break

        j1 += 1


    true_positive = true_count - j1

    while j2 < false_count:
        if false[j2] > i:
            j2 = j2 - 1
            break
        j2 += 1

    false_positive = false_count - j2
    sensitivity = float(true_positive)/true_count
    FPR = float(false_positive)/false_count
    sen_list.append(sensitivity)
    FPR_list.append(FPR)
plt.figure()
plt.plot(FPR_list, sen_list)
plt.savefig(args.fig)

"""
plt.figure()
bins = np.linspace(0, args.max_k, args.bin_size)
print tested_pair
plt.hist(true, bins, alpha=0.5, label='true', color = "r",normed=True)
plt.hist(false, bins, alpha=0.5, label='false', color = "b", normed=True)
plt.savefig(args.fig)
"""