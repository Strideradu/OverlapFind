import BlasrParse
import pickle
from Bio import SeqIO

blasr_large = {}
blasr_medium = {}
blasr_small = {}
score_threshold = 100


with open("D:/Data/20170312/sensitivety_large_overlap.m4") as f1:
    for line in f1:
        blasr_record = BlasrParse.BlasrParse(line)
        if blasr_record.score < score_threshold:
            pair = (blasr_record.target_id,blasr_record.id)
            if blasr_large.get(pair, False) is False:
                blasr_large[pair] = True

with open("D:/Data/20170312/sensitivety_medium_overlap.m4") as f1:
    for line in f1:
        blasr_record = BlasrParse.BlasrParse(line)
        if blasr_record.score < score_threshold:
            pair = (blasr_record.target_id,blasr_record.id)
            if blasr_medium.get(pair, False) is False:
                blasr_medium[pair] = True

with open("D:/Data/20170312/sensitivety_small_overlap.m4") as f1:
    for line in f1:
        blasr_record = BlasrParse.BlasrParse(line)
        if blasr_record.score < score_threshold:
            pair = (blasr_record.target_id,blasr_record.id)
            if blasr_small.get(pair, False) is False:
                blasr_small[pair] = True

def load_obj(filename ):
    with open(filename, 'rb') as f:
        return pickle.load(f)

num_test = 500
overlap_dict = load_obj("D:/Data/20170309/overlap.pkl")
fastq=SeqIO.index("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq")

unfoudn_pair = []
# print blasr_large
with open("D:/Data/20170309/blasr_overlap_not_found_score_0_1500pair.txt","w")as f:
    num_found = 0
    true_align = 0

    for pair in blasr_large:
        num_found += 1
        query = pair[0]
        target = pair[1]
        if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
            true_align += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

    print float(true_align) / num_found
    print >> f, ""

    num_found = 0
    true_align = 0

    for pair in blasr_medium:
        num_found += 1
        query = pair[0]
        target = pair[1]
        if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
            true_align += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

    print float(true_align) / num_found
    print >> f, ""

    num_found = 0
    true_align = 0

    for pair in blasr_small:
        num_found += 1
        query = pair[0]
        target = pair[1]
        if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
            true_align += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

    print float(true_align) / num_found
    print >> f, ""

