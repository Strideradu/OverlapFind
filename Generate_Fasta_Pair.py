from Bio import SeqIO

generate_num = 100
file_prefix = "D:/Data/20170706/FP_dustboth/fp_"
file_suffix = ".fasta"

query_fasta = "D:/Data/20170523/query_all.fasta"
target_fasta = "D:/Data/20170523/target_masked.fasta"
query_dict = SeqIO.index(query_fasta, "fasta")
target_dict = SeqIO.index(target_fasta, "fasta")
num = 0

with open("D:/Data/20170706/dustboth_fp_list.txt") as f:
    for line in f:
        num += 1
        pair = line.strip().split("\t")
        query_record = query_dict[pair[0]]
        target_record = target_dict[pair[1]]
        file_index = str(num).zfill(3)
        SeqIO.write(query_record, file_prefix +  "query_"+file_index + file_suffix, "fasta")
        SeqIO.write(target_record, file_prefix +  "target_"+file_index + file_suffix, "fasta")

        if num > generate_num:
            break

