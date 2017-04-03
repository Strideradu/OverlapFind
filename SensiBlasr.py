import BlasrParse
import pickle
from Bio import SeqIO

blasr_large = {}
blasr_medium = {}
blasr_small = {}
score_threshold = 1000

with open("D:/Data/20170402_BLASR_bestn/sensitivety_large_overlap_bestn_200.m4") as f1:
    negative_large = 0
    large_all_count = {}
    for line in f1:

        blasr_record = BlasrParse.BlasrParse(line)
        pair = (blasr_record.target_id, blasr_record.id)
        if large_all_count.get(pair, False) is False:
            negative_large += 1
            large_all_count[pair] = True

        if blasr_record.score < score_threshold:

            if blasr_large.get(pair, False) is False:
                blasr_large[pair] = True

print negative_large

with open("D:/Data/20170402_BLASR_bestn/sensitivety_medium_overlap_bestn_125.m4") as f1:
    negative_medium = 0
    medium_all_count = {}
    for line in f1:

        blasr_record = BlasrParse.BlasrParse(line)
        pair = (blasr_record.target_id, blasr_record.id)
        if medium_all_count.get(pair, False) is False:
            negative_medium += 1
            medium_all_count[pair] = True
        if blasr_record.score < score_threshold:

            if blasr_medium.get(pair, False) is False:
                blasr_medium[pair] = True

print negative_medium

with open("D:/Data/20170402_BLASR_bestn/sensitivety_small_overlap_bestn_125.m4") as f1:
    negative_small = 0
    small_all_count = {}
    for line in f1:

        blasr_record = BlasrParse.BlasrParse(line)
        pair = (blasr_record.target_id, blasr_record.id)
        if small_all_count.get(pair, False) is False:
            negative_small += 1
            small_all_count[pair] = True
        if blasr_record.score < score_threshold:
            pair = (blasr_record.target_id,blasr_record.id)
            if blasr_small.get(pair, False) is False:
                blasr_small[pair] = True

print negative_small

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
        if pair[0] != pair[1]:
            num_found += 1
            query = pair[0]
            target = pair[1]


            if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
                true_align += 1

            else:
                print >> f, pair[0] + "\t" + pair[1]
    print num_found
    # print "accuracy", float(true_align) / num_found
    print "sensitivity", float(true_align) / 610
    # print "FPR", (num_found - true_align)/float(negative_large - 610)
    print "FPR", (num_found - true_align) / float(42497- 610)
    print >> f, ""

    num_found = 0
    true_align = 0


    for pair in blasr_medium:
        if pair[0] != pair[1]:
            num_found += 1
            query = pair[0]
            target = pair[1]

            if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
                true_align += 1

            else:
                print >> f, pair[0] + "\t" + pair[1]
    print num_found
    # print "accuracy", float(true_align) / num_found
    print "sensitivity", float(true_align) / 588
    # print "FPR", (num_found - true_align) / float(negative_medium - 588)
    print "FPR", (num_found - true_align) / float(42499 - 588)
    print >> f, ""

    num_found = 0
    true_align = 0


    for pair in blasr_small:
        if pair[0] != pair[1]:
            num_found += 1
            query = pair[0]
            target = pair[1]

            if target in overlap_dict[query][0] or target in overlap_dict[query][1] or target in overlap_dict[query][2]:
                true_align += 1

            else:
                print >> f, pair[0] + "\t" + pair[1]
    print num_found
    # print "accuracy", float(true_align) / num_found
    print "sensitivity", float(true_align) /558
    # print "FPR", (num_found - true_align) / float(negative_small - 558)
    print "FPR", (num_found - true_align) / float(42498 - 558)
    print >> f, ""


