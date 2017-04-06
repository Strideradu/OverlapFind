import pickle
from Bio import SeqIO
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess
import argparse
import sys

"""
Script to test our method's sensitivity for discover true overlap pairs
The ground truth is divided into three categories as large, medium, small
we will test each categories and to sees the sensitivith
"""
parser = argparse.ArgumentParser()
parser.add_argument("accuracy", help="accuracy of the reads", type=float)
parser.add_argument("gap", help="gap rate", type=float)
parser.add_argument("rechain", help="length of kmer", type=int)
parser.add_argument("span", help="length of kmer", type=int)


try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

num_test = 500
#overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")
#fastq=SeqIO.index("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq")
#masked_fasta=SeqIO.index("D:/Data/20170116/filtered_15X_masked.fasta", "fasta")
overlap_dict = load_obj("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/overlap.pkl")
fastq=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/filtered_subreads_15X.fastq", "fastq")
masked_fasta=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170324_dust_group_hit/filtered_15X_masked.fasta", "fasta")

tested_query = {}
query = []
large_test = []
medium_test = []
small_test = []

in_large = {}
in_medium = {}
in_small = {}
true_pair = {}


for pacbio_id in overlap_dict:
    query.append(pacbio_id)
    large_overlap = overlap_dict[pacbio_id][0]
    medium_overlap = overlap_dict[pacbio_id][1]
    small_overlap = overlap_dict[pacbio_id][2]
    for overlap_seq in large_overlap:
        if tested_query.get(overlap_seq, False) is False:
            true_pair[(pacbio_id, overlap_seq)] = True
            if len(large_test)< num_test:
                if in_large.get(overlap_seq, False) is False:
                    large_test.append(overlap_seq)
                    in_large[overlap_seq] = True

    for overlap_seq in medium_overlap :
        if tested_query.get(overlap_seq, False) is False:
            true_pair[(pacbio_id, overlap_seq)] = True
            if len(medium_test)< num_test:
                if in_medium.get(overlap_seq, False) is False:
                    medium_test.append(overlap_seq)
                    in_medium[overlap_seq] = True

    for overlap_seq in small_overlap :
        if tested_query.get(overlap_seq, False) is False:
            true_pair[(pacbio_id, overlap_seq)] = True
            if len(small_test)< num_test :
                if in_small.get(overlap_seq, False) is False:
                    small_test.append(overlap_seq)
                    in_small[overlap_seq] = True

    tested_query[pacbio_id] = True
    # print pacbio_id
    if len(large_test)>=num_test and len(medium_test)>=num_test and len(small_test)>=num_test:
        break
#print len(query)


qual_seqs = {}
for test in [large_test, medium_test, small_test]:
    # num_pair = 0
    align_found = 0
    true_align = 0
    align_truth = 0
    tested = 0
    for query_seq in query:
        # print  query_seq
        large_overlap = overlap_dict[pacbio_id][0]
        medium_overlap = overlap_dict[pacbio_id][1]
        small_overlap = overlap_dict[pacbio_id][2]

        try:
            seq1 = qual_seqs[query_seq]
        except KeyError:
            record1 = fastq[query_seq]
            seq1 = QualitySeq(record1)
            qual_seqs[query_seq] = seq1

        for target_seq in test:
            # num_pair += 1
            if query_seq != target_seq:
                tested += 1

                try:
                    seq2 = qual_seqs[target_seq]
                except KeyError:
                    record2 = masked_fasta[target_seq]
                    seq2 = QualitySeq(record2)
                    qual_seqs[target_seq] = seq2

                process = DiagProcess(seq1, seq2)
                process.diag_points(9)
                chians = process.diag_chain(args.accuracy, args.gap)
                process.optimal_rechain(args.gap, args.rechain, args.span)
                if true_pair.get((query_seq, target_seq), False) is True or true_pair.get((target_seq, query_seq), False) is True:
                    align_truth += 1
                    if process.aligned:
                        align_found += 1
                        true_align += 1
                    """
                    else:
                        print >> f, query_seq + "\t" + target_seq
                    """
                else:
                    if process.aligned:
                        align_found += 1

                        # false positive align
                        print "False Positive" + query_seq + "\t" + target_seq
                        #sys.stdout.flush()

                print query_seq + "\t" + target_seq
                sys.stdout.flush()



    print "sensitivity: ", float(true_align) / align_truth
    print "FPR ", float(align_found - true_align) / (tested - align_truth)
    print "Found alignment: ", align_found
    print "Found True Align: ", true_align
    print "Ground Truth:", align_truth
    print "How many pair", tested
