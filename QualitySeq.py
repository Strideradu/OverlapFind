import numpy as np
from Kmer import Kmer

class QualitySeq(object):
    def __init__(self, fastq_record, insertion_record = None, deletion_record = None, substitution_record = None, del_tag = None, sub_tag = None):
        self.id = fastq_record.id
        self.seq = fastq_record.seq
        self.length = len(self.seq)
        self.kmer_count = {}

        Q = np.array(fastq_record.letter_annotations["phred_quality"])
        P = 10.0 ** (-Q / 10.0)
        self.quality = P

        if insertion_record:
            insertion_Q = np.array(insertion_record.letter_annotations["phred_quality"])
            insertion_P = 10.0 ** (-insertion_Q / 10.0)
            self.insertion_quality = insertion_P
        else:
            self.insertion_quality = None

        if deletion_record:
            deletion_Q = np.array(deletion_record.letter_annotations["phred_quality"])
            deletion_P = 10.0 ** (-deletion_Q / 10.0)
            self.deletion_quality = deletion_P
        else:
            self.deletion_quality = None

        if substitution_record:
            substitution_Q = np.array(substitution_record.letter_annotations["phred_quality"])
            substitution_P = 10.0 ** (-substitution_Q / 10.0)
            self.substitution_quality = substitution_P
        else:
            self.substitution_quality = None

        if del_tag:
            self.deletion_tag = del_tag.seq
        else:
            self.deletion_tag = None

        if sub_tag:
            self.subsitution_tag = sub_tag.seq
        else:
            self.subsitution_tag = None

    def generate_kmer(self, k, freq_dict, freq_threshold):
        kmers = []
        for i in range(self.length -k + 1):
            kmer = Kmer(str(self.seq[i:i+k]), freq_dict)
            if kmer.check_threshold(freq_threshold):
                kmers.append(kmer)

        return kmers

    def generate_good_kmer(self, k, freq_dict, accuracy_threshold, freq_threshold):
        qual_kmer  = []
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

                    accuracy = 1 - self.quality[i + n]
                    kmer_i += self.seq[n + i]
                    score += accuracy
                    k_i += 1

                    if len(kmer_i) == k:
                        # save the kmer, however we may add to the list later
                        kmer = Kmer(kmer_i, freq_dict)
                        if not (insertion_case):
                            if good_quality or kmer.check_threshold(freq_threshold):
                                qual_kmer.append(kmer)

                            break

                    if insertion_case:
                        if skipnext:
                            skipnext = False
                        else:
                            insert_kmer_i += self.seq[n + i]
                            insert_score += accuracy

                        if len(insert_kmer_i) == k:
                            insert_kmer = Kmer(insert_kmer_i, freq_dict)
                            qual_kmer.append(insert_kmer)

                            if good_quality or kmer.check_threshold(freq_threshold):
                                qual_kmer.append(kmer)

                            break

                    # print accuracy
                    if accuracy < accuracy_threshold:
                        good_quality = False

                        if not (insert_kmer_i):
                            insert_kmer_i = kmer_i
                            insert_score = score
                            insertion_case = True

                        if self.insertion_quality[i + n] > accuracy:
                            skipnext = True

                        else:
                            insertion_case = False
                            if kmer:
                                if good_quality or kmer.check_threshold(freq_threshold):
                                    qual_kmer.append(kmer)

                                break

                    n += 1

        return qual_kmer


