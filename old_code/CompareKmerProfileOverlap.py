# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 21:15:23 2017

@author: Nan
"""

import pickle
from Bio import SeqIO
import random
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import math


def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)
        
def count_score_kmer(record, k, quality_handle = None):
    kmer_list = []
    score_list = []
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
        score = np.sum(np.log(P[x:x + k]) - math.log(0.85))

        kmer_list.append(kmer_seq)
        score_list.append(score)
    return kmer_list, score_list
        
fastq=SeqIO.index("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170116_overlap_test/filtered_subreads_15X.fastq", "fastq")
id_list = list(SeqIO.parse("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170116_overlap_test/filtered_subreads_15X.fastq", "fastq"))
        
overlap_csv = "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/overlap_1.csv"

kmer_freq = load_obj("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/filtered_subreads_15X_kmer_freq.pkl")

abnormal_out = "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/abnormal_overlap.out"
    
true_freq = []
false_freq = []

true_score = []
false_score = []

true_overlap_count = []
false_overlap_count = []

forward_read_kmer_count = {}
reverse_read_kmer_count = {}
# print time.clock()
    
with open(overlap_csv) as f:
    with open(abnormal_out, "w") as f1: 
        for line in f:
            line = line.strip().split(",")
            seq_id = line[0]
            large = line[1].split()
            medium = line[2].split()
            small = line[3].split()
            overlap_list = large + medium + small
            overlap_num = len(large + medium)
            
            seq = fastq[seq_id]
            query_kmer, query_score = count_score_kmer(seq, 11)
            # query_count = Counter(query_kmer)
            
            n = 0
            picked_id = []
            
            while n < overlap_num:
                # select non overlap sequence and generate false sample
                sample = random.choice(id_list)
                if sample.id in picked_id or sample.id in overlap_list:
                    continue
                else:
                    forward_sample_dict = forward_read_kmer_count.get(sample.id,False)
                    reverse_sample_dict = reverse_read_kmer_count.get(sample.id, False)
                    if not(forward_sample_dict):
                        target_kmer, target_score = count_score_kmer(sample, 11)
                        forward_sample_dict = Counter(target_kmer)
                        forward_read_kmer_count[sample.id] = forward_sample_dict
                        target_kmer, target_score = count_score_kmer(sample.reverse_complement(), 11)
                        reverse_sample_dict = Counter(target_kmer)
                        reverse_read_kmer_count[sample.id] = reverse_sample_dict
                    
                    read_forward_overlap = []
                    read_forward_score = []
                    read_reverse_overlap = []
                    read_reverse_score = []
                    
                    forward_overlap_count = 0
                    for i, kmer in enumerate(query_kmer):
                        if forward_sample_dict.get(kmer, False):
                            read_forward_overlap.append(kmer_freq[kmer])
                            read_forward_score.append(query_score[i])
                            forward_overlap_count += 1
                            
                    reverse_overlap_count = 0
                    for i, kmer in enumerate(query_kmer):
                        if forward_sample_dict.get(kmer, False):
                            read_reverse_overlap = [].append(kmer_freq[kmer])
                            read_reverse_score .append(query_score[i])
                            reverse_overlap_count += 1
                        
                    if forward_overlap_count >= reverse_overlap_count:
                        false_freq.extend(read_forward_overlap)
                        false_score.extend(read_forward_score)
                        false_overlap_count.append(forward_overlap_count)
                        
                    else:
                        false_freq.extend(read_reverse_overlap)
                        false_score.extend(read_reverse_score)
                        false_overlap_count.append(reverse_overlap_count)
                    #print overlap_count
                    picked_id.append(sample.id)
                    n += 1
        
               
            for target_id in large + medium:

                forward_sample_dict = forward_read_kmer_count.get(target_id,False)
                reverse_sample_dict = reverse_read_kmer_count.get(target_id, False)
                if not(forward_sample_dict):
                    tagrte_seq = fastq[target_id]
                    target_kmer, target_score = count_score_kmer(tagrte_seq, 11)
                    forward_sample_dict = Counter(target_kmer)
                    forward_read_kmer_count[target_id] = forward_sample_dict
                    target_kmer, target_score = count_score_kmer(tagrte_seq.reverse_complement(), 11)
                    reverse_sample_dict = Counter(target_kmer)
                    reverse_read_kmer_count[target_id] = reverse_sample_dict
            
              
                read_forward_overlap = []
                read_forward_score = []
                read_reverse_overlap = []
                read_reverse_score = []
                
                forward_overlap_count = 0
                for i, kmer in enumerate(query_kmer):
                    if forward_sample_dict.get(kmer, False):
                        read_forward_overlap.append(kmer_freq[kmer])
                        read_forward_score.append(query_score[i])
                        forward_overlap_count += 1
                        
                reverse_overlap_count = 0
                for i, kmer in enumerate(query_kmer):
                    if forward_sample_dict.get(kmer, False):
                        read_reverse_overlap = [].append(kmer_freq[kmer])
                        read_reverse_score .append(query_score[i])
                        reverse_overlap_count += 1
                    
                if forward_overlap_count >= reverse_overlap_count:
                    true_freq.extend(read_forward_overlap)
                    true_score.extend(read_forward_score)
                    true_overlap_count.append(forward_overlap_count)
                    overlap_count = forward_overlap_count
                    
                else:
                    true_freq.extend(read_reverse_overlap)
                    true_score.extend(read_reverse_score)
                    true_overlap_count.append(reverse_overlap_count)
                    overlap_count = reverse_overlap_count
                
                # detect really bad overlap case
                if overlap_count * 11 <= 0.25 * len(seq.seq):
                    print >> f1, ",".join([seq_id, target_id, str(overlap_count)])
                
            #print true_overlap_count, false_overlap_count
        
# print time.clock()

plt.figure()
bins = np.linspace(-50, 5, 100)
plt.hist(true_score, bins, alpha=0.5, label='true', color = "r")
plt.hist(false_score, bins, alpha=0.5, label='false', color = "b")
plt.savefig("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/11mer_score.png")

plt.figure()
bins = np.linspace(0, 100, 120)
plt.hist(true_freq, bins, alpha=0.5, label='true', color = "r")
plt.hist(false_freq, bins, alpha=0.5, label='false', color = "b")
plt.savefig("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/11mer_spectrum.png")                
            
plt.figure()
bins = np.linspace(0, 1000, 200)
plt.hist(true_overlap_count, bins, alpha=0.5, label='true', color = "r")
plt.hist(false_overlap_count, bins, alpha=0.5, label='false', color = "b")
plt.savefig("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170120_overlapkmer/11mer_overlap_count_for_read.png")           
        
        
    