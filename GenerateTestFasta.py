from Bio import SeqIO
import random

ratio = float(4)/70
print ratio
num_query_job = 12

query_id = [[], [], [], [], [], [], [], [], [], [], [], []]
print query_id
records =  list(SeqIO.parse("D:/Data/20170116/filtered_subreads_15X.fastq","fastq"))
masked_index = SeqIO.index("D:/Data/20170116/filtered_15X_masked.fasta","fasta")
num_read = 0

print len(records)
for i in range(len(records)):
    r = random.random()
    if r < ratio:
        num_read += 1
        query_id[num_read % num_query_job].append(i)

normal_read = []
masked_read = []
file_dir ="D:/Data/20170523/query_"
file_ext = ".fasta"
for i in range(num_query_job):
    output = []
    file_path = file_dir + str(i).zfill(2) + file_ext
    for j in range(len(query_id[i])):
        record = records[query_id[i][j]]
        output.append(record)
        normal_read.append(record)
        masked_read.append(masked_index[record.id])

    SeqIO.write(output, file_path, "fasta")

SeqIO.write(normal_read, "D:/Data/20170523/query_all.fasta", "fasta")
SeqIO.write(masked_read, "D:/Data/20170523/target_masked.fasta", "fasta")
