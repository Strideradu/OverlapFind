# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 12:13:17 2017

@author: Nan
"""

import os
import sys
from Bio import SeqIO
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import Counter
import quality_kmer


def count_score_kmer(record, k, kmer_list, quality_handle = None):
    if quality_handle == None:
        fastq_record = record
    else:
        fastq_record = quality_handle[record.id]
    Q = np.array(fastq_record.letter_annotations["phred_quality"])
    P = 10.0 ** (-Q / 10.0)
    # P is the error probability

    # total_kmers = len(record.seq) - k + 1
    # first assemble dict of kmer counts

    for x in range(len(record.seq) + 1 - k):
        kmer_seq = record.seq[x:x + k]
        score = np.sum(P[x:x + k])

        kmer_list.append(kmer_seq)
"""
fastq = sys.argv[1]
if os.path.exists(fastq):
    print os.path.basename(fastq)
""" 
fastq = "D:/Data/20170116/filtered_subreads_15X.fastq"   
records = SeqIO.parse(fastq, "fastq")
insertion_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_insertion_first1k.fastq", "fastq")
deletion_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_deletion_first1k.fastq", "fastq")
substitution_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_substitution_first1k.fastq", "fastq")
deletion_tag_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_deletion_tag_first1k.fasta", "fasta")
substitution_tag_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_substitution_tag_first1k.fasta", "fasta")


"""
output = sys.argv[2]
if os.path.exists(output):
    print os.path.basename(output)
""" 
output = "D:/Data/20170116/filtered_subreads_15X_kmer.fasta" 

kmer_list = []

for record in records:
    insertion_rec = insertion_dict[record.id]
    deletion_rec = deletion_dict[record.id]
    substitution_rec = substitution_dict[record.id]
    sub_tag_rec = substitution_tag_dict[record.id]
    del_tag_rec = deletion_tag_dict[record.id]

    qual_record = quality_kmer.QualitySeq(record, insertion_rec, deletion_rec, substitution_rec, del_tag_rec, sub_tag_rec)
    record_kmer = qual_record.generate_quality_kmer(11)    
    
    
kmer_num = Counter(record_kmer)

kmer_fasta = [] 

for kmer in kmer_num:
    kmer_fasta.append(SeqRecord(Seq(kmer),id=kmer + "_" + str(kmer_num[kmer]), description=""))
    

SeqIO.write(kmer_fasta, output, "fasta")