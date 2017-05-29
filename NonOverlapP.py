from ProbFunc import r_randommatch_probability
import numpy as np

p = r_randommatch_probability(50, 9, 8000, 8000)
np.save("/mnt/home/dunan/Job/2016/201605_align_noisy_long-reads/20170529_large_p_matrix/p_matrix.npy", p)
