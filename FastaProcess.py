import os
from Bio import SeqIO
from quality_seq import QualitySeq

class FastaProcess(object):
    def __init__(self, file_path, type = "fasta"):
        self.file = file_path
        self.type = type

    def quality_kmer_fasta(self, output_path, insertion_path, deletion_path = None, substitution_path = None, deletion_tag_path = None, substitution_tag_path = None):
        """

        :param output_path: path to output the fasta file for kmer
        :param insertion_path: insertion quality fastq file
        :param deletion_path:
        :param substitution_path:
        :param deletion_tag_path:
        :param substitution_tag_path:
        :return:
        """
        if self.type == "fastq":
            records = SeqIO.parse(self.file, self.type)
            insertion_dict = SeqIO.index(insertion_path)

            if deletion_path:
                deletion_dict = SeqIO.index(deletion_path)
            else:
                deletion_dict = {}

            if substitution_path:
                substitution_dict = SeqIO.index(substitution_path)
            else:
                substitution_dict = {}

            if deletion_tag_path:
                deletion_tag_dict = SeqIO.index(deletion_tag_path)
            else:
                deletion_tag_dict = {}

            if substitution_tag_path:
                substitution_tag_dict = SeqIO.index(substitution_tag_path)
            else:
                substitution_tag_dict = {}

            for record in records:
                insertion_rec = insertion_dict[record.id]
                deletion_rec = deletion_dict.get(record.id, None)
                substitution_rec = substitution_dict.get(record.id, None)
                sub_tag_rec = substitution_tag_dict.get(record.id, None)
                del_tag_rec = deletion_tag_dict.get(record.id, None)

