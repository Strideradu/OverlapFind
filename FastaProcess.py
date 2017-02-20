import os
from Bio import SeqIO
from QualitySeq import QualitySeq
import pickle
import argparse
import sys


class FastaProcess(object):
    def __init__(self, file_path, type = "fasta"):
        self.file = file_path
        self.type = type
        self.kmer_freq = None

    def load_freq(self, dict_path):
        with open(dict_path, 'rb') as f:
            self.kmer_freq = pickle.load(f)

    def kmer_fasta(self, k, output_path, freq_threshold=0):
        """

        :param k:  size of kmer
        :param output_path:
        :param freq_threshold:
        :return:
        """
        with open(output_path, "w") as output_handle:
            records = SeqIO.parse(self.file, self.type)
            for record in records:
                qual_record = QualitySeq(record)
                qual_kmer = qual_record.generate_kmer(k, self.kmer_freq, freq_threshold)

                for kmer in qual_kmer:
                    kmer.to_fasta(output_handle, "fasta")


    def quality_kmer_fasta(self, k, output_path, insertion_path, deletion_path = None, substitution_path = None, deletion_tag_path = None, substitution_tag_path = None, accuracy_threshold = 0.75, freq_threshold = 0):
        """

        :param output_path: path to output the fasta file for kmer
        :param insertion_path: insertion quality fastq file
        :param deletion_path:
        :param substitution_path:
        :param deletion_tag_path:
        :param substitution_tag_path:
        :return:
        """
        with open(output_path, "w") as output_handle:
            if self.type == "fastq":
                records = SeqIO.parse(self.file, self.type)
                insertion_dict = SeqIO.index(insertion_path, "fastq")

                if deletion_path:
                    deletion_dict = SeqIO.index(deletion_path, "fastq")
                else:
                    deletion_dict = {}

                if substitution_path:
                    substitution_dict = SeqIO.index(substitution_path, "fastq")
                else:
                    substitution_dict = {}

                if deletion_tag_path:
                    deletion_tag_dict = SeqIO.index(deletion_tag_path, "fastq")
                else:
                    deletion_tag_dict = {}

                if substitution_tag_path:
                    substitution_tag_dict = SeqIO.index(substitution_tag_path, "fastq")
                else:
                    substitution_tag_dict = {}

                for record in records:
                    insertion_rec = insertion_dict[record.id]
                    deletion_rec = deletion_dict.get(record.id, None)
                    substitution_rec = substitution_dict.get(record.id, None)
                    sub_tag_rec = substitution_tag_dict.get(record.id, None)
                    del_tag_rec = deletion_tag_dict.get(record.id, None)

                    qual_record = QualitySeq(record, insertion_rec, deletion_rec, substitution_rec, del_tag_rec, sub_tag_rec)
                    qual_kmer = qual_record.generate_good_kmer(k, self.kmer_freq, accuracy_threshold, freq_threshold)

                    for kmer in qual_kmer:
                        kmer.to_fasta(output_handle, "fasta")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("k", help="length of kmer", type=int)
    parser.add_argument("template_fastq", help="the fastq file that contain sequence we want to extract kmer", type=str,  nargs=1)
    parser.add_argument("-insert_Q", dest="insertion_fastq", help="the fastq file that contain insertion quality", type=str,  nargs=1)
    parser.add_argument("output_fasta", help="the fasta file to store the extracted kmer", type=str,  nargs=1)
    parser.add_argument("-freq_dict", dest="dict", help="path that store frequency dict", type=str,  nargs=1)
    parser.add_argument("-freq", dest="freq", help="the threshold to filter frequency", type=int,  nargs=1)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    if os.path.exists(args.template_fastq[0]):
        print os.path.basename(args.template_fastq[0])

    process = FastaProcess(args.template_fastq[0], "fastq")

    if args.dict:
        process.load_freq(args.dict[0])

    if args.insertion_fastq:

        if args.freq:
            process.quality_kmer_fasta(args.k, args.output_fasta[0], args.insertion_fastq[0], freq_threshold=args.freq[0])
        else:
            process.quality_kmer_fasta(args.k, args.output_fasta[0], args.insertion_fastq[0])

    else:
        if args.freq:
            process.kmer_fasta(args.k,  args.output_fasta[0], freq_threshold=args.freq[0])
        else:
            process.kmer_fasta(args.k, args.output_fasta[0])

if __name__ == '__main__':
    main()