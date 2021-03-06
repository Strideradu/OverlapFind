import argparse
from Bio import SeqIO
import sys

sys.path.append("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/OverlapFind/")
import PseudoBloomFilter as bf
import QuerySeqPro as QuerySeq
import ProbFunc

parser = argparse.ArgumentParser()
parser.add_argument("query", help="path of query file", type=str)
parser.add_argument("target", help="path of target file", type=str)
parser.add_argument("k", help="kmer", type=int)
parser.add_argument("accuracy", help="accuracy of the reads", type=float)
parser.add_argument("chain_size", help="length of kmer for a cluster to be considered as true align", type=int)
parser.add_argument("span_coefficient", help="must be non-zero float", type=float)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

query_seq = list(SeqIO.parse(args.query, "fasta"))
target_seq = list(SeqIO.parse(args.target, "fasta"))
print "# kmer size:", args.k
print "# match rate:", args.accuracy
print "# chian size threshold:", args.chain_size


L = ProbFunc.statistical_bound_of_waiting_time(args.accuracy, args.k)
print "# Group Distance:", L

for record2 in target_seq:

    for record1 in query_seq:
        filter = bf.PseudoBloomFilter(record2, args.k, L)
        filter.generate_filter()

        query = QuerySeq.QuerySeq(record1)
        query.check_kmer(filter)
        query.cluster_hits(size_threshold=args.chain_size, debug=False, group_hit=args.span_coefficient)

        if query.aligned:
            print record1.id + "\t" + record2.id + "\t" + "True"

