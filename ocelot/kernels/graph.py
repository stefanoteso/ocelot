# -*- coding: utf8 -*-

import numpy as np

from .base import _Kernel

class DiffusionKernel(_Kernel):
    """The diffusion kernel between graph nodes [Kondor02]."""
    def __init__(self, adjmatrix, beta = 1.0, **kwargs):
        n, m = adjmatrix.shape
        assert n == m
        self._adjmatrix = adjmatrix
        self._beta = beta
        super(DiffusionKernel, self).__init__(range(n), **kwargs)
    def _compute_all(self):
        import scipy.linalg as la
        # XXX this does not work for weighted adjacency matrices
        a = self._adjmatrix
        d = np.identity(a.shape[0]) * np.sum(a, axis = 1)
        e = la.expm(self._beta * (a - d))

        # being approximate, matrix exponentiation can give slightly asymmetric
        # results; here we fix them if they are small enough.
        asymmetry = np.linalg.norm(e.T - e, ord = "fro")
        if asymmetry > 0.0:
            assert asymmetry < 1e-10, "the gods are playful today"
            e = 0.5 * (e.T + e)
        return e
