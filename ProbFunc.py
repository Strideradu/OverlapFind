""""
functions for probability from YASS
"""
import numpy as np

def statistical_bound_of_waiting_time(p, k, alpha = 0.05):
    """

    :param p:  average accuracy of read
    :param k:
    :param alpha: statistical bound
    :return:statistical bound of waiting time
    """
    p_k = float(p)**(k)
    qp_k = (1-p)*p_k
    x  = 0
    sum_0_xk1 = 0.00
    sum = 0

    last_k_prob = np.zeros((k+1), dtype = np.double)

    while sum < 1 - alpha:
        if x < k:
            last_k_prob[x % (k + 1)] = 0.00
        elif x == k:
            last_k_prob[x % (k + 1)] = p_k
        else:
            sum_0_xk1 += last_k_prob[x % (k + 1)]
            last_k_prob[x % (k + 1)] = qp_k * (1 - sum_0_xk1)

        sum += last_k_prob[x % (k + 1)]
        x += 1

    return x

def randomwalk_probability_of_pos(pI, L):
    """

    :param pI: rates for insertion and deletions
    :param L: bound of waiting time
    :return:
    """
    P = L
    a = pI * 0.50
    b = 1 - pI

    u = np.zeros((2 * L +1), dtype = np.double)
    f = np.zeros((2 * L +1), dtype = np.double)
    t = np.zeros((2 * L +1), dtype = np.double)


    for i in range(2 * L + 1):
        u[i] = 0
    for i in range(2 * L + 1):
        f[i] = 0

    f[0] = 1.00
    u[0] = a
    u[1] = b
    u[2] = a

    while P > 0:
        if P & 1:
            for i in range(2 * L + 1):
                t[i] = 0
                for j in range(i + 1):
                    t[i] += u[i-j] * f[j]

            s = t
            t = f
            f = s

        for i in range(2 * L + 1):
            t[i] = 0
            for j in range(i + 1):
                t[i] += u[i - j] * u[j]

        s = t
        t = u
        u = s

        P>>=1

    return f

def statistical_bound_of_randomwalk(pI, L, alpha = 0.05):
    """

    :param pI: rates for insertion and deletions
    :param L: bound of waiting time
    :param alpha: statistical bound
    :return:
    """
    RDW_Bound = randomwalk_probability_of_pos(pI, L)
    sum = RDW_Bound[L]
    bound = 0

    while sum < (1-alpha) and bound < L:
        bound += 1
        sum += RDW_Bound[L - bound] + RDW_Bound[L + bound]

    return bound

def r_matches_probability(r, k, p, L):
    """

    :param r: How many kmer hits?
    :param k: size of kmer
    :param p: accuracy
    :param L: sequence leangth
    :return: probability
    """
    # DP problem, need to fill a matrix of L*r
    p_matrix = np.zeros((L + 1, r + 1), dtype=np.double)

    # fill the [x, r = 0] of the matrix
    p_k = float(p) ** (k)
    qp_k = (1 - p) * p_k
    sum_0_xk1 = 0 # store the x-k-1 term of sun
    last_k_prob = np.zeros((k + 1), dtype=np.double)
    sum = 0

    for x in range(L + 1 - k):


        if x < k:
            last_k_prob[x % (k + 1)] = 0.00
        elif x == k:
            last_k_prob[x % (k + 1)] = p_k
        else:
            sum_0_xk1 += last_k_prob[x % (k + 1)]
            last_k_prob[x % (k + 1)] = qp_k * (1 - sum_0_xk1)
            # last_k_prob[x % (k + 1)] = p_k * (1 - sum_0_xk1)

        sum += last_k_prob[x % (k + 1)]

        if x < k:
            p_matrix[x][0] = 1.0

        p_matrix[x + k][0] = 1 - sum

    # special case x = k , r = 1 p = p**k


    for i in range(L + 1):
        for j in range(1, r+1):
            if i > k:
                # if i < k, the p is always 0
                p_matrix[i][j] = p_matrix[i-1][j] + qp_k*(p_matrix[i-k-1][j-1] - p_matrix[i-k-1][j])

    return p_matrix[L][r]


def r_randommatch_probability(r, size_k, L1, L2, p = 0.25):
    """

    :param r: number of kmer hit
    :param k: k of kmer
    :param L1: length of read1
    :param L2: length of read2
    :param p: random match rata for 1 base
    :return:
    """
    p_matrix = np.zeros((L1 + 1, L2 + 1, r + 1), dtype=np.double)

    # initialize p_matrix
    p_k = float(p) ** (size_k)
    qp_k = (1 - p) * p_k
    sum_0_xk1 = 0  # store the x-k-1 term of sun
    last_k_prob = np.zeros((size_k + 1), dtype=np.double)
    sum = 0

    for x in range(L1 + 1):
        if x < size_k:
            last_k_prob[x % (size_k + 1)] = 0.00
        elif x == size_k:
            last_k_prob[x % (size_k + 1)] = p_k
        else:
            sum_0_xk1 += last_k_prob[x % (size_k + 1)]
            last_k_prob[x % (size_k + 1)] = qp_k * (1 - sum_0_xk1)
            # last_k_prob[x % (k + 1)] = p_k * (1 - sum_0_xk1)

        sum += last_k_prob[x % (size_k + 1)]

        p_matrix[x][size_k][0] = 1 - sum

    ##print p_matrix[k][k][0]
    #print p_matrix[11][size_k][0]

    sum_0_xk1 = 0  # store the x-k-1 term of sun
    last_k_prob = np.zeros((size_k + 1), dtype=np.double)
    sum = 0

    for x in range(L2 + 1):
        if x < size_k:
            last_k_prob[x % (size_k + 1)] = 0.00
        elif x == size_k:
            last_k_prob[x % (size_k + 1)] = p_k
        else:
            sum_0_xk1 += last_k_prob[x % (size_k + 1)]
            last_k_prob[x % (size_k + 1)] = qp_k * (1 - sum_0_xk1)
            # last_k_prob[x % (k + 1)] = p_k * (1 - sum_0_xk1)

        sum += last_k_prob[x % (size_k + 1)]

        p_matrix[size_k][x][0] = 1 - sum
    #print p_matrix[size_k][11][0]
    #print p_matrix[k][k][0]
    for i in range(L1 + 1):
        for j in range(L2 + 1):
            if i > size_k and j > size_k:
                sum_diag = 0
                for k in range(1, size_k):
                    sum_diag += (p_matrix[i-k-1][j-k-1][0])*p**(k)*(1-p)
                p_matrix[i][j][0] = (p_matrix[i - 1][j][0] - p_matrix[i-1][j-1][0]) * (1 - p) \
                                    + (p_matrix[i][j-1][0] - p_matrix[i-1][j-1][0]) * (1 - p) \
                                    + p_matrix[i - 1][j - 1][0] * (1 - p) \
                                    + sum_diag
                if p_matrix[i][j][0] < 0:
                    p_matrix[i][j][0] = 0.0

    print "test",p_matrix[L1][L2][0]
    #print p_matrix[L1 -30][L2 - 30][0]

    for k in range(1, r + 1):
        for i in range(1, L1 + 1):
            for j in range(1, L2 + 1):

                if i > size_k  and j > size_k:

                    p_matrix[i][j][k] = (p_matrix[i - 1][j][k] - p_matrix[i-1][j-1][k]) \
                                        + (p_matrix[i][j-1][k] - p_matrix[i-1][j-1][k]) \
                                        + p_matrix[i-1][j-1][k] \
                                        + qp_k * (p_matrix[i - size_k - 1][j - size_k - 1][k - 1])

                else:
                    p_matrix[i][j][k] = (p_matrix[i - 1][j][k] - p_matrix[i - 1][j - 1][k])  \
                                        + (p_matrix[i][j - 1][k] - p_matrix[i - 1][j - 1][k]) \
                                        + p_matrix[i - 1][j - 1][k]
    #print p_matrix[L1/2][L2/2][r]
    #print p_matrix[L1 - 50][L2 - 50][1]
    for k in range(r):
        print k, p_matrix[L1][L2][k]

    return p_matrix[L1][L2][r]

if __name__ == '__main__':


    #L = statistical_bound_of_waiting_time(0.8, 9)
    #print L
    #print statistical_bound_of_randomwalk(0.15, L)

    print r_matches_probability(30, 13, 0.85*0.85, 5000)
    # so if we change accuracy to 1- accuracy, we can calculate the probability of have r matches given two non-overlap reads?
    # if we have p(1-0.85*0.85), then we have r matches given two non overlap reads is ((1/3)**k*(x-k))**r*p + (1/4)**(kr)*L

    #print r_randommatch_probability(50, 9, 2000, 5000, p=0.25)
    #print r_randommatch_probability(5, 9, 5000, 2000, p=0.25)
    print r_randommatch_probability(5, 9, 50, 100, p = 0.25)
    print r_randommatch_probability(5, 9, 100, 50, p = 0.25)


