from __future__ import division

import numpy as np

def random_k_sparse(n, k, sigma2):
    '''
    Generate a k sparse signal of length n whose entries are
    drawn at random from a zero-mean Normal distribution with
    variance sigma2.
    '''

    supp = np.random.choice(n, size=k, replace=False)
    val = np.random.randn(k) * np.sqrt(sigma2)

    x = np.zeros(n)
    x[supp] = val

    return x, val, supp
