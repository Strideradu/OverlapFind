import pickle
from Bio import SeqIO
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess

"""
Script to test our method's sensitivity for discover true overlap pairs
The ground truth is divided into three categories as large, medium, small
we will test each categories and to sees the sensitivith
"""

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

num_test = 500
overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")
fastq=SeqIO.index("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq")

tested = {}
query = []
large_test = []
medium_test = []
small_test = []

in_large = {}
in_medium = {}
in_small = {}

for pacbio_id in overlap_dict:
    query.append(pacbio_id)
    large_overlap = overlap_dict[pacbio_id][0]
    medium_overlap = overlap_dict[pacbio_id][1]
    small_overlap = overlap_dict[pacbio_id][2]
    for overlap_seq in large_overlap:
        if len(large_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_large.get(overlap_seq, False) is False:
                large_test.append(overlap_seq)
                in_large[overlap_seq] = True

    for overlap_seq in medium_overlap :
        if len(medium_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_medium.get(overlap_seq, False) is False:
                medium_test.append(overlap_seq)
                in_medium[overlap_seq] = True

    for overlap_seq in small_overlap :
        if len(small_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_small.get(overlap_seq, False) is False:
                small_test.append(overlap_seq)
                in_small[overlap_seq] = True

    tested[pacbio_id] = True

    if len(large_test)>=num_test and len(medium_test)>=num_test and len(small_test)>=num_test:
        break

with open("D:/Data/20170309/overlap_not_found_1500pair.txt","w")as f:
    for test in [large_test, medium_test, small_test]:
        # num_pair = 0
        align_found = 0
        true_align = 0
        align_truth = 0
        for query_seq in query:
            large_overlap = overlap_dict[pacbio_id][0]
            medium_overlap = overlap_dict[pacbio_id][1]
            small_overlap = overlap_dict[pacbio_id][2]
            for target_seq in test:
                # num_pair += 1
                record1 = fastq[query_seq]
                record2 = fastq[target_seq]
                seq1 = QualitySeq(record1)
                seq2 = QualitySeq(record2)
                process = DiagProcess(seq1, seq2)
                process.diag_points(9)
                chians = process.diag_chain(0.75, 0.2)
                process.rechain(0.2)
                if target_seq in large_overlap or target_seq in medium_overlap or target_seq in small_overlap:
                    align_truth += 1
                    if process.aligned:
                        align_found += 1
                        true_align += 1
                    else:
                        print >> f, query_seq + "\t" + target_seq
                else:
                    if process.aligned:
                        align_found += 1

        print "sensitivity: ", float(true_align) / align_truth
        print "accuracy: ", float(true_align) / align_found
        print align_truth
        print >> f, ""