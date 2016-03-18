# -*- coding: utf-8 -*-

from os.path import join, basename
import numpy as np
import pandas as pd
from itertools import product
from collections import defaultdict, Counter
from glob import glob

from ocelot.kernels import *
from ocelot.microarray import read_pcl, GeneExprKernel
from ocelot.services import InterProTSV

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

class InterProKernel(Kernel):
    """A simple domain kernel built around InterPro.

    Parameters
    ----------
    ps : collection
        Ordered list of protein IDs.
    path : str
        Path to the directory holding the InterPro files.
    mode : str, optional
        One of "match", "count", "evalue", defaults to "match".
    allowed_sources : collection or None, optional
        List of allowed domain providers, defaults to all of them.
    default_score : float
        Score to use when no E-value is provided, defaults to 1.0.
    """
    def __init__(self, ps, path, mode="match", allowed_sources=None,
                 default_score=1.0, *args, **kwargs):
        if not mode in ("match", "count", "evalue"):
            raise ValueError("invalid mode '{}'".format(mode))
        self._path = path
        self._mode = mode
        self._allowed_sources = allowed_sources
        self._default_score = default_score
        super(InterProKernel, self).__init__(ps, *args, **kwargs)

    def _to_score(self, evalue):
        if evalue is None or evalue <= 0.0:
            return self._default_score
        return -np.log(evalue)

    def _compute_all(self):
        parser = InterProTSV()

        all_hits, num_missing = [], 0
        for p in self._entities:
            try:
                path = join(self._path, "{}.tsv.txt".format(p))
                domain_to_evalue = parser.read(path, self._allowed_sources)
            except IOError, e:
                domain_to_evalue = {}
                num_missing += 1

            if self._mode == "match":
                hits = set(domain_to_evalue.keys())
            elif self._mode == "count":
                hits = dict(Counter(domain_to_evalue.keys()))
            elif self._mode == "evalue":
                hits = {domain: self._to_score(evalue)
                        for domain, evalue in domain_to_evalue.iteritems()}
            all_hits.append(hits)

        if num_missing > 0:
            print "no interpro domains for '{}/{}' proteins" \
                    .format(num_missing, len(self))

        if self._mode == "match":
            return SetKernel(all_hits).compute()
        else:
            return SparseLinearKernel(all_hits).compute()
