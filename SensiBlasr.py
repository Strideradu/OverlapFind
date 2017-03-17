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

tested = {}
large_test = []
medium_test = []
small_test = []
for pacbio_id in overlap_dict:
    large_overlap = overlap_dict[pacbio_id][0]
    medium_overlap = overlap_dict[pacbio_id][1]
    small_overlap = overlap_dict[pacbio_id][2]
    for overlap_seq in large_overlap:
        if len(large_test)< num_test and tested.get(overlap_seq, False) is False:
            large_test.append((pacbio_id, overlap_seq))

    for overlap_seq in medium_overlap :
        if len(medium_test)< num_test and tested.get(overlap_seq, False) is False:
            medium_test.append((pacbio_id, overlap_seq))

    for overlap_seq in small_overlap :
        if len(small_test)< num_test and tested.get(overlap_seq, False) is False:
            small_test.append((pacbio_id, overlap_seq))

    tested[pacbio_id] = True

    if len(large_test)>=num_test and len(medium_test)>=num_test and len(small_test)>=num_test:
        break


unfoudn_pair = []
# print blasr_large
with open("D:/Data/20170309/blasr_overlap_not_found_score_0_1500pair.txt","w")as f:
    num_pair = 0
    num_found = 0
    for pair in large_test:
        if blasr_large.get(pair, False) is True:
            num_found += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

        num_pair += 1
    print float(num_found) / num_pair
    print >> f, ""

    num_pair = 0
    num_found = 0
    for pair in medium_test:
        if blasr_medium.get(pair, False) is True:
            num_found += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

        num_pair += 1
    print float(num_found) / num_pair
    print >> f, ""

    num_pair = 0
    num_found = 0
    for pair in small_test:
        if blasr_small.get(pair, False) is True:
            num_found += 1

        else:
            print >> f, pair[0] + "\t" + pair[1]

        num_pair += 1
    print float(num_found) / num_pair
    print >> f, ""


