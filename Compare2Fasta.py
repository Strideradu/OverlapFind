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
parser.add_argument("accuracy", help="accuracy of the reads", type=float)
parser.add_argument("gap", help="gap rate", type=float)
parser.add_argument("rechain", help="length of kmer", type=int)
parser.add_argument("span_coefficient", help="must be non-zero float", type=float)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

query_seq = list(SeqIO.parse(args.query, "fasta"))
target_seq = list(SeqIO.parse(args.target, "fasta"))

L = ProbFunc.statistical_bound_of_waiting_time(args.accuracy, args.k)
delta = ProbFunc.statistical_bound_of_randomwalk(args.gap, L)

for record1 in query_seq:

    for record2 in target_seq:
        # print record1.id, record2.id
        seq1 = QualitySeq(record1)
        seq2 = QualitySeq(record2)

        process = DiagProcess(seq1, seq2)
        process.diag_points(args.k)
        process.diag_group_hit(L, delta)
        process.optimal_rechain(args.gap, args.rechain, args.span_coefficient)

        if process.aligned:
            print record1.id + "\t" + record2.id + "\t" + "True"