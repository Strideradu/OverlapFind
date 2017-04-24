import argparse
from Bio import SeqIO
import sys
from QualitySeq import QualitySeq
from DiagProcess import DiagProcess

parser = argparse.ArgumentParser()
parser.add_argument("query", help="path of query file", type=str)
parser.add_argument("target", help="path of target file", type=str)
parser.add_argument("k", help="kmer", type=int)
parser.add_argument("accuracy", help="accuracy of the reads", type=float)
parser.add_argument("gap", help="gap rate", type=float)
parser.add_argument("rechain", help="length of kmer", type=int)
parser.add_argument("span", help="length of kmer", type=int)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

query_seq = SeqIO.parse(args.query, "fasta")
target_seq = SeqIO.index(args.target, "fasta")
target_ids = target_seq.keys()

for record1 in query_seq:

    for record2_id in target_ids:
        print record1.id, record2_id
        seq1 = QualitySeq(record1)
        record2 = target_seq[record2_id]
        seq2 = QualitySeq(record2)

        process = DiagProcess(seq1, seq2)
        process.diag_points(args.k)
        process.diag_chain(args.accuracy, args.gap)
        process.optimal_rechain(args.gap, args.rechain, args.span)

        if process.aligned:
            print record1.id + "\t" + record2.id