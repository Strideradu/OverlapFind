# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:01:10 2017
Open m4 file, determine overlap of all the reads from alignments
@author: Nan
"""
import ParseOutput
from intervaltree import Interval, IntervalTree
import pickle

m4_align = "D:/Data/20170116/mapped.m4"
m4_overlap_output = "D:/Data/20170126/overlap.csv"

m4_records = {}
tree = IntervalTree()


def save_obj(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


with open(m4_align) as f:
    for line in f:
        m4_record = ParseOutput.BlasrRecord(line)
        aligned_len = m4_record.query_end - m4_record.query_start

        if aligned_len > 0.9 * m4_record.query_len:
            old_record = m4_records.get(m4_record.id, False)
            if not (old_record) or m4_record.score < old_record.score:
                m4_records[m4_record.id] = m4_record

m4_lists = m4_records.keys()
for m4_id in m4_lists:
    m4_record = m4_records[m4_id]
    tree[m4_record.target_start:m4_record.target_end] = m4_id

overlap_dict = {}
for m4_id in m4_lists:
    m4_record = m4_records[m4_id]
    length = m4_record.query_len
    large = []
    medium = []
    small = []
    overlap_list = list(tree.search(m4_record.target_start, m4_record.target_end))
    # print overlap_list
    for overlap_rec in overlap_list:
        if overlap_rec.data != m4_id:
            x = range(m4_record.target_start, m4_record.target_end)
            y = range(overlap_rec.begin, overlap_rec.end)
            # print x
            # print y
            # print set(x) & set(y)
            ovelap_len = len(set(x) & set(y))
            overlap_frac = float(ovelap_len) / length
            if overlap_frac >= 0.5:
                large.append(overlap_rec.data)
            elif overlap_frac >= 0.25:
                medium.append(overlap_rec.data)
            elif overlap_frac > 0:
                small.append(overlap_rec.data)

    overlap_dict[m4_id] = [large, medium, small]

save_obj(overlap_dict, "D:/Data/20170309/overlap.pkl")
