# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 19:50:38 2017
Generate frequency dict
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
import pickle

def save_obj(obj, filename ):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

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
        
fastq = "D:/Data/20170116/filtered_subreads_15X.fastq"   
records = SeqIO.parse(fastq, "fastq")

kmer_list = []

for record in records:
    record_kmer = count_score_kmer(record, 11, kmer_list)
    
kmer_num = Counter(kmer_list)

save_obj(kmer_num,"D:/Data/20170116/filtered_subreads_15X_kmer_freq.pkl")