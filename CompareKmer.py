import argparse
from Bio import SeqIO
import sys
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess
import ProbFunc

parser = argparse.ArgumentParser()
parser.add_argument("query", help="path of query file", type=str)
parser.add_argument("target", help="path of target file", type=str)
parser.add_argument("k", help="kmer", type=int)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

query_seq = list(SeqIO.parse(args.query, "fasta"))
target_seq = list(SeqIO.parse(args.target, "fasta"))

for record1 in query_seq:

    for record2 in target_seq:
        # print record1.id, record2.id
        seq1 = QualitySeq(record1)
        seq2 = QualitySeq(record2)

        process = DiagProcess(seq1, seq2)
        process.diag_points(args.k)
        fw_num = len(process.fw_points)
        rc_num = len(process.rc_points)
        hits_num = max(fw_num, rc_num)

        if hits_num != 0:
            print record1.id + "\t" + record2.id + "\t" + str(hits_num) + "\t" + str(len(record1.seq)) + "\t" + str(len(record2.seq))