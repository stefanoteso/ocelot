# -*- coding: utf-8 -*-

from sklearn.utils import check_random_state

def permute(l, rng=None):
    rng = check_random_state(rng)

    perm = list(rng.permutation(len(l)))
    return [l[perm[i]] for i in range(len(l))]
