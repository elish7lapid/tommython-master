import numpy as np


def estimate_time(N1, N2, alpha=1):
    """
    Assumptions: sequence1 and sequence 2 are of te same length,
    over the alphabet {a,c,g,t,-}
    We use a Jukes-Cantor matrix of the following form:
        -4alpha alpha alpha alpha alpha
        alpha -4alpha alpha alpha alpha
        alpha alpha -4alpha alpha alpha
        alpha alpha alpha -4alpha alpha
        alpha alpha alpha alpha -4alpha
    Our Jukes-Cantor model has a 5X5 matrix, because it compares between sequences with gaps.
    Our model assumes that each position is changing independently, according to a continuous time Markov chain
    With the above Jukes-Cantor Rate Matrix.
    :param N1: will be the number of equal places in the alignment
    :param N2: will be the number of distinct places
    :param alpha: the Jukes-Cantor parameter
    :return: the maximum likelihood estimation of the time passed between the 2 sequences.
    """
    # N1 will be the number of equal places in the alignment
    # N2 will be the number of distinct places

    if 4*N1 <= N2:
        return np.inf
    return (float(1)/float(5*alpha))*np.log((4*(N1+N2))/(4*N1-N2))


def probability(a, b, t, alpha=1):
    """
    the probability that a will become b after time t, under Jukes-Cantor 5X5 model
    :param a: the former letter
    :param b: the later letter
    :param t: time passed
    :param alpha: the Jukes Cantor parameter
    :return: the probability that a will become b after time t, under Jukes-Cantor 5X5 model
    """
    if t == np.inf:
        return 0.2
    if a == b:
        return 0.2*(1 + 4*np.exp(-5*alpha*t))
    else:
        return 0.2*(1 - np.exp(-5*alpha*t))
