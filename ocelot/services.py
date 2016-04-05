# -*- coding: utf-8 -*-

import os, errno, hashlib
import numpy as np
import multiprocessing as mp
import itertools as it
from ocelot.sequences import read_fasta, write_fasta, AMINOACIDS


def iterate_csv(path, num_skip = 0, **kwargs):
    """Thin wrapper around ``csv.DictReader`` with a couple more options."""
    import csv
    with open(path, "rU") as fp:
        for row_as_dict in csv.DictReader(fp, **kwargs):
            if num_skip > 0:
                num_skip -= 1
                continue
            yield row_as_dict
