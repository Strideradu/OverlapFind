from Bio import SeqIO
import pickle
from matplotlib import pyplot as plt


def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


overlap_dict = load_obj("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170317_ROC/overlap.pkl")



file_list = [
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/9mer/compare_kmer_all_9mer.out",
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/11mer/compare_kmer_all_11mer.out",
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/12mer/compare_kmer_all_12mer.out",
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/13mer/compare_kmer_all_13mer.out",
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/14mer/compare_kmer_all_14mer.out",
    "/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/15mer/compare_kmer_all_15mer.out"
]

plt.figure()
for file in file_list:
    true = []
    false = []
    with open(file)as result_file:
        for line in result_file:

            # check first symble of seq id
            if line[0] == "m":
                line_sp = line.strip().split()
                query_id = line_sp[0]
                target_id = line_sp[1]
                num_hits = int(line_sp[2])

                if query_id != target_id:
                    if overlap_dict.get(query_id, False) and (
                                        target_id in overlap_dict[query_id][0] or target_id in overlap_dict[query_id][
                                    1] or target_id in overlap_dict[query_id][2]):
                        true.append(num_hits)

                    else:
                        false.append(num_hits)

    true.sort()
    false.sort()
    true_count = len(true)
    false_count = len(false)
    max_value = max(true[-1], false[-1])

    sen_list = []
    FPR_list = []

    j1 = 0
    j2 = 0
    for i in range(max_value):
        while j1 < true_count:
            if true[j1] > i:
                j1 = j1 - 1
                if j1 < 0:
                    j1 = 0
                break

            j1 += 1

        true_positive = true_count - j1

        while j2 < false_count:
            if false[j2] > i:
                j2 = j2 - 1
                if j2 < 0:
                    j2 = 0
                break
            j2 += 1

        false_positive = false_count - j2
        sensitivity = float(true_positive) / true_count
        FPR = float(false_positive) / 2*(347361 - 342)
        sen_list.append(sensitivity)
        FPR_list.append(FPR)

    plt.plot(FPR_list, sen_list, label = file.split("_")[-1])
plt.xlim(xmin=-0.0005, xmax=0.1005)
plt.savefig("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170821_kmer_retest/kmer_compare.png")
