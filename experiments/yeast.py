# -*- coding: utf-8 -*-

import os
import numpy as np
import itertools as it
from collections import Counter
from glob import glob

from ocelot.kernels import *
from ocelot.services import _cls, iterate_csv, PCL, InterProTSV

class SGDGeneExpressionKernel(Kernel):
    """A yeast-specific correlation kernel on gene expression data.

    It uses the microarray files located at::

        $src/yeast/microarray/*.pcl

    The data can be obtained from:

        `<ftp://downloads.yeastgenome.org/expression/microarray/all_spell_exprconn_pcl.tar.gz>`_

    :param p_to_i: map between ORF feature names to indices in the kernel matrix.
    :param src: path to the data source directory.
    :param experiments: list of microarray experiments to use.
    """
    def __init__(self, p_to_i, src, experiments, *args, **kwargs):
        self._wc = os.path.join(src, "SGD", "microarray", "*", "*.pcl")
        self._experiments = experiments
        super(SGDGeneExpressionKernel, self).__init__(p_to_i, *args, **kwargs)

    def _get_pcl_paths(self):
        for path in glob(self._wc):
            if not any(exp in path for exp in self._experiments):
                continue
            yield path

    def _compute_all(self):
        pcl = PCL()
        p_to_i = self._entities
        self._matrix = np.zeros((len(self), len(self)), dtype=self.dtype)
        for path in self._get_pcl_paths():
            print _cls(self), ": processing '{}'".format(path)
            p_to_levels, num_conditions = pcl.read(path)
            exp_levels, num_missing = [], 0
            for p, i in sorted(p_to_i.iteritems(), key = lambda p_i: p_i[1]):
                try:
                    p_levels = p_to_levels[p].astype(self.dtype)
                except:
                    p_levels = np.zeros((num_conditions,), dtype=self.dtype)
                    num_missing += 1
                exp_levels.append(p_levels)
            if num_missing > 0:
                 print _cls(self), ": '{}' has no measurements for '{}/{}' proteins" \
                                    .format(path, num_missing, len(self))
            self._matrix += CorrelationKernel(exp_levels, do_normalize=False).compute()
        return self._matrix

class YeastProteinComplexKernel(Kernel):
    """A yeast-specific diffusion kernel on protein complexes.

    It uses the protein complex file located at::

        $src/yeast/ppi/CYC2008_complex.tab

    :param p_to_i: map between ORF feature names to indices in the kernel matrix.
    :param src: path to the data source directory.
    :param gamma: parameter of the diffusion kernel.
    """
    def __init__(self, p_to_i, src, *args, **kwargs):
        self._path = os.path.join(src, "yeast", "ppi", "CYC2008_complex.tab")
        super(YeastProteinComplexKernel, self).__init__(p_to_i, *args, **kwargs)

    def _read_complex_to_orf(self):
        FIELDS = ("ORF", "_", "COMPLEX")
        complex_to_orfs = {}
        for row in iterate_csv(self._path, delimiter = "\t", num_skip = 1,
                               fieldnames = FIELDS):
            orf, cplex = row["ORF"], row["COMPLEX"]
            if not cplex in complex_to_orfs:
                complex_to_orfs[cplex] = set()
            complex_to_orfs[cplex].add(orf)
        return complex_to_orfs

    def _compute_all(self):
        p_to_i = self._entities

        complex_to_orfs = self._read_complex_to_orf()

        adjmatrix = np.zeros((len(self), len(self)))
        for _, orfs in complex_to_orfs.iteritems():
            orf_indices = [p_to_i[orf] for orf in orfs if orf in p_to_i]
            for i, j in it.product(orf_indices, orf_indices):
                if i != j:
                    adjmatrix[i,j] = 1
        return DiffusionKernel(adjmatrix).compute()

class InterProKernel(Kernel):
    """A simple domain kernel built around InterPro.

    :param ps: list of protein IDs.
    :param cache_path: path to the directory holding the interpro output.
    :param allowed_sources: list of allowed domain providers (default: ``None``).
    :param mode: one of ``"match"``, ``"count"``, ``"evalue"`` (default: ``"match"``).
    :param default_score: default score to use when no evalue is provided (default: ``1.0``).
    """
    def __init__(self, ps, cache_path, allowed_sources = None,
                 use_evalue = False, default_score = 1.0, *args, **kwargs):
        if not mode in ("match", "count", "evalue"):
            raise ValueError("invalid mode '{}'".format(mode))
        self._cache_path = cache_path
        self._allowed_sources = allowed_sources
        self._mode = mode
        self._default_score = default_score
        super(InterProKernel, self).__init__(ps, *args, **kwargs)

    def _path(self, p):
        return os.path.join(self._cache_path, "interpro",
                            "{}.tsv.txt".format(p))

    def _compute_iprscan_hits(self):
        # WRITEME
        pass

    def _to_score(self, evalue):
        if evalue is None or evalue <= 0.0:
            return self._default_score
        return -np.log(evalue)

    def _compute_all(self):
        self._compute_iprscan_hits()
        parser = InterProTSV()

        all_hits, num_missing = [], 0
        for p in self._entities:
            try:
                domain_to_evalue = parser.read(self._path(p),
                                               self._allowed_sources)
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
            print _cls(self), ": no interpro domains for '{}/{}' proteins" \
                                .format(num_missing, len(self))

        if self._mode == "match":
            return SetKernel(all_hits).compute()
        else:
            return SparseLinearKernel(all_hits).compute()
