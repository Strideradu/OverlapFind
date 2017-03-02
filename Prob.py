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
    sum_0_xk1 = 0.00;
    sum = 0;

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

    return x;

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
        u[i] = 0;
    for i in range(2 * L + 1):
        f[i] = 0;

    f[0] = 1.00;
    u[0] = a;
    u[1] = b;
    u[2] = a;

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

if __name__ == '__main__':
     L = statistical_bound_of_waiting_time(0.85, 6)
     print L
     print statistical_bound_of_randomwalk(0.12, L)


