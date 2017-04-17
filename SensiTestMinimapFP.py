import pickle
from Bio import SeqIO
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("k", help="kmer size", type=int)
parser.add_argument("accuracy", help="accuracy of the reads", type=float)
parser.add_argument("gap", help="gap rate", type=float)
parser.add_argument("rechain", help="length of kmer", type=int)
parser.add_argument("span", help="length of kmer", type=int)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

print args.k, args.accuracy, args.gap, args.rechain, args.span

fastq=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/filtered_subreads_15X.fastq", "fastq")
masked_fasta=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170324_dust_group_hit/filtered_15X_masked.fasta", "fasta")

qual_seqs = {}

fp_num = 0
fp_found = 0
with open("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170414_test_minimap_missing_case/minimap_15X_k13w5_fp_align.out") as f:
    for line in f:
        line = line.strip()
        pair = line.split("\t")
        query_id = pair[0]
        target_id = pair[1]

        try:
            seq1 = qual_seqs[query_id]
        except KeyError:
            record1 = fastq[query_id]
            seq1 = QualitySeq(record1)
            qual_seqs[query_id] = seq1

        try:
            seq2 = qual_seqs[target_id]
        except KeyError:
            record2 = masked_fasta[target_id]
            seq2 = QualitySeq(record2)
            qual_seqs[target_id] = seq2

        fp_num += 1
        process = DiagProcess(seq1, seq2)
        process.diag_points(args.k)
        process.diag_chain(args.accuracy, args.gap)
        process.optimal_rechain(args.gap, args.rechain, args.span)
        if process.aligned:
            fp_found += 1

        print query_id + "\t" + target_id + "\t finished"
        sys.stdout.flush()

print "found aligned pair in missing dataset", fp_num - fp_found
print "minimpa missed hit", fp_num
sys.stdout.flush()