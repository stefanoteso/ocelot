# -*- coding: utf-8 -*-

from os.path import join, basename
import numpy as np
import pandas as pd
from itertools import product
from collections import defaultdict, Counter
from glob import glob

from ocelot.kernels import *
from ocelot.microarray import read_pcl, GeneExprKernel

class SGDGeneExprKernel(GeneExprKernel):
    """A yeast-specific gene expression kernel.

    It uses the PCL files located at::

        ${src}/SGD/microarray/${experiments}/*.pcl

    The data can be obtained from [1]_.

    Parameters
    ----------
    gene_to_i : dict
        Map from gene names to indices in the Gram matrix. There should be
        a key for each target gene.
    src : str
        WRITEME
    experiments : list
        WRITEME
    """
    def __init__(self, gene_to_i, src, experiments, *args, **kwargs):
        microarrays = self._load(src, experiments)
        super(SGDGeneExprKernel, self).__init__(gene_to_i, microarrays,
                                                *args, **kwargs)

    def _load(self, src, experiments):
        microarrays = []
        basepath = join(src, "SGD", "microarray")
        for path in glob(join(basepath, "*")):
            experiment = basename(path)
            if not experiment in experiments:
                continue
            for pcl in glob(join(path, "*.pcl")):
                print "loading '{}'".format(pcl)
                microarrays.append(read_pcl(pcl))
        return microarrays

class YeastProteinComplexKernel(Kernel):
    """A yeast-specific diffusion kernel on protein complexes.

    It uses the protein complex file located at::

        ${src}/yeast/ppi/CYC2008_complex.tab

    Parameters
    ----------
    p_to_i : dict
        Map ORF feature name -> index.
    src : str
        Path to the data directory.
    beta : float
        Diffusion parameter.
    """
    def __init__(self, p_to_i, src, beta=1.0, *args, **kwargs):
        self._path = join(src, "yeast", "ppi", "CYC2008_complex.tab")
        self._beta = beta
        super(YeastProteinComplexKernel, self).__init__(p_to_i, *args, **kwargs)

    def _compute_all(self):
        df = pd.read_csv(self._path, sep="\t")
        complex_to_orfs = defaultdict(set)
        for _, row in df.iterrows():
            complex_to_orfs[row.Complex].add(row.ORF)

        p_to_i = self._entities
        matrix = np.zeros((len(self), len(self)), dtype=self.dtype)
        for orfs_in_complex in complex_to_orfs.itervalues():
            indices = [p_to_i[orf] for orf in orfs_in_complex if orf in p_to_i]
            matrix[np.ix_(indices, indices)] = 1

        return DiffusionKernel(matrix, beta=self._beta).compute()
