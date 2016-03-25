# -*- coding: utf-8 -*-

import subprocess
import numpy as np
from os import makedirs
from os.path import isdir, abspath
from sklearn.utils import check_random_state

def checkdir(path):
    return path and isdir(abspath(path))

def quietmkdir(path):
    try:
        makedirs(cache)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise RuntimeError("can not create cache directory '{}': {}" \
                                .format(cache, e))

def validate(valid, values):
    if values is None:
        values = valid
    else:
        for value in values:
            if not value in valid:
                raise ValueError("invalid value '{}'".format(value))
    return values

def permute(l, rng=None):
    rng = check_random_state(rng)

    perm = list(rng.permutation(len(l)))
    return [l[perm[i]] for i in range(len(l))]

def ispsd(matrix, tol=1e-6):
    return (np.abs(matrix.T - matrix) < tol).all() and \
           np.all(np.linalg.eigvalsh(matrix) >= 0.0)

def split_tr_ts(array, indices0, indices1=None):
    if indices1 is None:
        indices1 = sorted(set(range(len(array.shape[0]))) - set(indices0))
    if array.ndim == 1:
        return array[indices0], array[indices1]
    elif array.ndim == 2:
        assert array.shape[0] == array.shape[1]
        return array[np.ix_(indices1, indices1)], \
               array[np.ix_(indices0, indices1)]
    else:
        raise ValueError("invalid ndim")

def run_binary(path, args, shell=True):
    pipe = subprocess.Popen(path + " " + " ".join(args), shell=shell,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = pipe.communicate()
    return pipe.wait(), out, err
