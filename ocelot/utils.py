# -*- coding: utf-8 -*-

import numpy as np
from sklearn.utils import check_random_state

def permute(l, rng=None):
    rng = check_random_state(rng)

    perm = list(rng.permutation(len(l)))
    return [l[perm[i]] for i in range(len(l))]

def split_tr_ts(array, indices0, indices1):
    if array.ndim == 1:
        return array[indices0], array[indices1]
    elif array.ndim == 2:
        return array[np.ix_(tr_indices, tr_indices)], \
               array[np.ix_(ts_indices, tr_indices)]
    else:
        raise ValueError("invalid ndim")
