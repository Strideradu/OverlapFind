import pickle
from Bio import SeqIO

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
in_query = {}
large_test = []
in_large = {}
medium_test = []
in_medium = {}
small_test = []
in_small = {}
for pacbio_id in overlap_dict:
    if in_query.get(pacbio_id, False) is False:
        query.append(fastq[pacbio_id])
        in_query[pacbio_id] = True
    large_overlap = overlap_dict[pacbio_id][0]
    medium_overlap = overlap_dict[pacbio_id][1]
    small_overlap = overlap_dict[pacbio_id][2]
    for overlap_seq in large_overlap:
        if len(large_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_large.get(overlap_seq, False) is False:
                large_test.append(fastq[overlap_seq])
                in_large[overlap_seq] = True

    for overlap_seq in medium_overlap :
        if len(medium_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_medium.get(overlap_seq, False) is False:
                medium_test.append(fastq[overlap_seq])
                in_medium[overlap_seq] = True

    for overlap_seq in small_overlap :
        if len(small_test)< num_test and tested.get(overlap_seq, False) is False:
            if in_small.get(overlap_seq, False) is False:
                small_test.append(fastq[overlap_seq])
                in_small[overlap_seq] = True

    tested[pacbio_id] = True

    if len(large_test)>=num_test and len(medium_test)>=num_test and len(small_test)>=num_test:
        break

with open("D:/Data/20170312/sensitivety_query.fasta","w") as query_fasta:
    SeqIO.write(query, query_fasta, "fasta")
with open("D:/Data/20170312/sensitivety_large_overlap.fasta","w") as large_overlap_fasta:
    SeqIO.write(large_test, large_overlap_fasta , "fasta")
with open("D:/Data/20170312/sensitivety_medium_overlap.fasta","w") as medium_overlap_fasta:
    SeqIO.write(medium_test, medium_overlap_fasta , "fasta")
with open("D:/Data/20170312/sensitivety_small_overlap.fasta","w") as small_overlap_fasta:
    SeqIO.write(small_test, small_overlap_fasta , "fasta")