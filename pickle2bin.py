#!/usr/bin/env python3

import sys
import pickle
import numpy as np

if len(sys.argv) != 3:
    print("Usage: {} <pickle input> <binary output>".format(sys.argv[0]))
    sys.exit(1)

print("loading...")
with open(sys.argv[1], "rb") as fp:
    # XXX the encoding thing works around python 2/3 incompatibility, see:
    #
    #   http://stackoverflow.com/questions/11305790/pickle-incompatability-of-numpy-arrays-between-python-2-and-3
    #
    try:
        matrix = pickle.load(fp)
    except:
        matrix = pickle.load(fp, encoding='latin1')

print("storing upper triangle...")
with open(sys.argv[2], "wb") as fp:
    info = np.array([matrix.shape[0]], dtype=np.int32)
    info.tofile(fp)
    for i, row in enumerate(matrix):
        row[i:].astype(np.float32).tofile(fp)
