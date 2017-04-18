from Bio import SeqIO

fastq=SeqIO.index("D:/Data/20170116/filtered_subreads_15X.fastq", "fastq")
masked_fasta=SeqIO.index("D:/Data/20170116/filtered_15X_masked.fasta", "fasta")

total_len = 0
num = 0
search_space = 0
long_read = 0
short_read = 0
with open("D:/Data/20170414/minimap_15X_k13w5_missing_align.out") as f:
    for line in f:
        line = line.strip()
        pair = line.split("\t")
        query_id = pair[0]
        target_id = pair[1]
        record1 = fastq[query_id]
        record2 = masked_fasta[target_id]
        num += 2
        len_1 = len(record1)
        len_2 = len(record2)
        total_len +=  len_1 + len_2
        search_space += len_1*len_2
        if len_1 > len_2:
            long_read += len_1
            short_read += len_2
        else:
            long_read += len_2
            short_read += len_1

print "average len: " + str(float(total_len)/num)
print "average long len: " + str(2.0*long_read/num)
print "average short len: " + str(2.0*short_read/num)
print "average search space: " + str(2.0*search_space/num)

total_len = 0
num = 0
search_space = 0
long_read = 0
short_read = 0
with open("D:/Data/20170414/minimap_15X_k13w5_fp_align.out") as f:
    for line in f:
        line = line.strip()
        pair = line.split("\t")
        query_id = pair[0]
        target_id = pair[1]
        record1 = fastq[query_id]
        record2 = masked_fasta[target_id]
        num += 2
        len_1 = len(record1)
        len_2 = len(record2)
        total_len += len_1 + len_2
        search_space += len_1 * len_2
        if len_1 > len_2:
            long_read += len_1
            short_read += len_2
        else:
            long_read += len_2
            short_read += len_1

print "average len: " + str(float(total_len)/num)
print "average long len: " + str(2.0*long_read/num)
print "average short len: " + str(2.0*short_read/num)
print "average search space: " + str(2.0*search_space/num)