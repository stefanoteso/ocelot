#!/usr/bin/env python2

import sys
import cPickle as pickle

if len(sys.argv) != 3:
    print "Usage: {} <input pickle file> <output txt file>".format(sys.argv[0])
    sys.exit(1)

print "loading..."
with open(sys.argv[1], "rb") as fp:
    matrix = pickle.load(fp)

print "writing..."
with open(sys.argv[2], "wb") as fp:
    fp.write("{} {}\n".format(matrix.shape[0], matrix.shape[1]))
    for row in matrix:
        fp.write(" ".join("{0:.4f}".format(elem) for elem in row) + "\n")

print "done"
