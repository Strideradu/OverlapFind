from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import numpy as np
import random

random.seed(0)
accuracy_threshold = 0.75
insertion_threshold = 0.3
frequency_threshold = 12

class QualitySeq(object):
    def __init__(self, fastq_record, insertion_record, deletion_record, substitution_record, del_tag, sub_tag):
        self.seq = fastq_record.seq
        self.length = len(self.seq)
        self.kmer_count = {}

        Q = np.array(fastq_record.letter_annotations["phred_quality"])
        P = 10.0 ** (-Q / 10.0)
        self.quality = P

        insertion_Q = np.array(insertion_record.letter_annotations["phred_quality"])
        insertion_P = 10.0 ** (-insertion_Q / 10.0)
        self.insertion_quality = insertion_P

        deletion_Q = np.array(deletion_record.letter_annotations["phred_quality"])
        deletion_P = 10.0 ** (-deletion_Q / 10.0)
        self.deletion_quality = deletion_P

        substitution_Q = np.array(substitution_record.letter_annotations["phred_quality"])
        substitution_P = 10.0 ** (-substitution_Q / 10.0)
        self.substitution_quality = substitution_P

        self.deletion_tag = del_tag.seq
        self.subsitution_tag = sub_tag.seq

    def generate_quality_kmer(self, k, insertion_kmer = False):
        qual_kmer = []
        for i, base in enumerate(self.seq):
            # print i
            if i < self.length - k + 1:
                k_i = 0
                n = 0
                kmer_i = ""
                score = 0
                while k_i < k and i + n < self.length:

                    accuracy = 1 - self.quality[i+n]
                    # print accuracy
                    if accuracy > accuracy_threshold:
                        kmer_i += self.seq[n + i]
                        score += accuracy
                        k_i += 1
                        n += 1
                    else:
                        if insertion_kmer:
                            if self.insertion_quality[i + n] > accuracy:
                                # insertion case
                                n += 1
                                continue
                            else:
                                break
                        else:
                            break


                if len(kmer_i) == k:
                    qual_kmer.append(kmer_i)
        return qual_kmer

    def freq_check(kmer, freq_dict, kmer_list):
        if freq_dict.get(kmer, 0) > frequency_threshold:
            kmer_list.append(kmer)
        
    
    def generate_good_kmer(self, k, insertion_kmer = False, freq_dict):
        qual_kmer = []
        for i, base in enumerate(self.seq):
            # print i
            if i < self.length - k + 1:
                k_i = 0
                n = 0
                kmer_i = ""
                insert_kmer_i = None
                score = 0
                good_quality = True
                insertion_case = False     
                skipnext = False
                kmer = None
                
                while i + n < self.length:

                    accuracy = 1 - self.quality[i+n]
                    
                    kmer_i += self.seq[n + i]
                    score += accuracy
                    k_i += 1
                    
                    if len(kmer_i) == k:
                        kmer = kmer_i
                        if not(insertion_case):
                            if not(good_quality):
                                freq_check(kmer, freq_dict, qual_kmer)
                            else:
                                qual_kmer.append(kmer)
                                
                            break                    
                    
                    if insertion_case:
                        if skipnext:
                            skipnext = False
                        else: 
                            insert_kmer_i += self.seq[n + i]
                            insert_score += accuracy
                            
                        if len(insert_kmer_i) == k:
                            qual_kmer.append(insert_kmer_i)
                            
                            if not(good_quality):
                                freq_check(kmer, freq_dict, qual_kmer)
                            else:
                                qual_kmer.append(kmer)
                            
                            break
                            
                    
                            
                    
                    # print accuracy
                    if accuracy < accuracy_threshold:
                        good_quality = False
                        
                        if not(insert_kmer_i):
                            insert_kmer_i = kmer_i
                            insert_score = score
                            insertion_case = True
                            
                        if self.insertion_quality[i + n] > insertion_threshold:
                            skipnext = True
                            
                        else:
                            insertion_case = False
                            if kmer:
                                if not(good_quality):
                                    freq_check(kmer, freq_dict, qual_kmer)
                                else:
                                    qual_kmer.append(kmer)
                                
                            
                    n += 1    
                        
                    
        return qual_kmer        
        
        
    def generate_quality_kmer_2(self, k):
        qual_kmer = []
        for i, base in enumerate(self.seq):

            if i < self.length - k + 1:
                k_i = 0
                n = 0
                kmer_i = ""
                while k_i < 14:
                    accuracy = 1 - self.quality[i+n]
                    max_qual = accuracy + self.insertion_quality[i+n] + self.deletion_quality[i+n] + \
                               self.substitution_quality[i+n]

                    m = random.random()

                    # print accuracy

                    if 0 <= m <= accuracy:
                        kmer_i += self.seq[n + i]
                        k_i += 1

                    else:
                        m2 = random.random()
                        if 0 < m2 <= self.deletion_quality[i + n]:
                            if self.deletion_tag[i + n] != "N":
                                kmer_i += self.deletion_tag[i + n] + self.seq[n + i]
                                k_i += 2
                            else:
                                self.deletion_tag[i + n] + self.seq[n + i]
                                k_i += 1
                        max3 = self.insertion_quality[i + n] + self.substitution_quality[i + n]
                        m3 = max3 * random.random()
                        if 0 < m3 <= self.insertion_quality[i + n]:
                            continue
                        else:
                            kmer_i += self.subsitution_tag[i + n]
                            k_i += 1
                    n += 1

                qual_kmer.append(kmer_i)
        return qual_kmer


    def generate_no_error_seq(self):
        no_error_seq = ""
        for i, base in enumerate(self.seq):
            accuracy = 1 - self.quality[i]      # this is general accuracy
            m = random.random()

            # since quality < 1, so we can use power of quality to make error rate smaller
            del_prob = (self.deletion_quality[i])**1

            if 0 <= m <= del_prob:
                if self.deletion_tag[i] != "N":
                    no_error_seq += self.deletion_tag[i] + base
                else:
                    no_error_seq += base

            m2 = random.random()
            m3 = random.random()

            sub_prob = (self.substitution_quality[i])**1
            #sub_prob = 0
            in_prob = self.insertion_quality[i]

            if 0 <= m2 <= sub_prob :
                if (0 <= m3 <= in_prob and m2 > m3) or m3 > in_prob:
                    no_error_seq += self.subsitution_tag[i]

            if 0 <= m3 <= in_prob :
                if (0 <= m2 <= sub_prob and m3 > m2) or m2 > sub_prob:
                    continue

            if m > del_prob and m2 > sub_prob and m3 > in_prob:
                no_error_seq += base

        return no_error_seq



if __name__ == "__main__":
    # specify the sequence file and the other infromation file
    records = SeqIO.parse("D:/Data/20161201/filtered_subreads_first100.fastq", "fastq")
    insertion_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_insertion_first1k.fastq", "fastq")
    deletion_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_deletion_first1k.fastq", "fastq")
    substitution_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_substitution_first1k.fastq", "fastq")
    deletion_tag_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_deletion_tag_first1k.fasta", "fasta")
    substitution_tag_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_substitution_tag_first1k.fasta", "fasta")

    # output quality to csv file
    quality_output = open("D:/Data/20161129/seq2_quality.csv", "w")

    # for fasta output
    result = []
    j = 1
    for record in records:

        insertion_rec = insertion_dict[record.id]
        deletion_rec = deletion_dict[record.id]
        substitution_rec = substitution_dict[record.id]
        sub_tag_rec = substitution_tag_dict[record.id]
        del_tag_rec = deletion_tag_dict[record.id]

        qual_record = QualitySeq(record, insertion_rec, deletion_rec, substitution_rec, del_tag_rec, sub_tag_rec)

        seq_list = list(qual_record.seq)
        del_list = list(qual_record.deletion_tag)
        sub_list = list(qual_record.subsitution_tag)
        print >> quality_output, ",".join(seq_list)
        print >> quality_output, ",".join(map(str,1 - qual_record.quality))
        print >> quality_output, ",".join(map(str, qual_record.insertion_quality))
        print >> quality_output, ",".join(map(str, qual_record.deletion_quality))
        print >> quality_output, ",".join(del_list)
        print >> quality_output, ",".join(map(str, qual_record.substitution_quality))
        print >> quality_output, ",".join(sub_list)

        result.append( SeqRecord(Seq(qual_record.generate_no_error_seq(), generic_dna), id =record.id, description=""))
        j += 1
        # following are test code to generate kmer
        """
        kmer_original = qual_record.generate_quality_kmer(14)
        print len(kmer_original)
        kmer = list(set(kmer_original))
        print len(kmer)
        """
    SeqIO.write(result, "D:/Data/20161210/first_100_random_generate.fasta", "fasta")