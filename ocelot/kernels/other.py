# -*- coding: utf8 -*-

import numpy as np
from itertools import product

from ocelot.kernels import Kernel

class EmpiricalKernel(Kernel):
    """Empirical kernel map over an underlying kernel [1]_.

    Given a fixed set of patterns :math:`(z_1, \\ldots, z_m)` and an input
    pattern :math:`x`, the empirical kernel map is defined as:

    .. math::

        \\varphi(x) := (k(x,z_1), ..., k(x,z_m))

    for any given sub-kernel :math:`k`.

    Parameters
    ----------
    indices : ordered list
        Indices of anchor patterns in the subkernel.
    subkernel : np.ndarray of shape (n, n) or Kernel
        The subkernel.

    References
    ----------
    .. [1] Scholkopf et al., "Input Space Versus Feature Space in Kernel-Based
           Methods", 1999.
    """
    def __init__(self, indices, subkernel, **kwargs):
        super(EmpiricalKernel, self).__init__(indices, **kwargs)
        self._subkernel = subkernel

    def _compute_all(self):
        indices = self._entities

        try:
            submatrix = self._subkernel.compute()
        except AttributeError:
            submatrix = self._subkernel

        phi = submatrix[indices].astype(self.dtype)
        return np.dot(phi, phi.T)

class PairwiseKernel(Kernel):
    """The pairwise kernel, defined as:

    .. math::
        k[(i,j),(n,m)] := k[i,n]k[j,m] + k[i,m]k[j,n]

    or:

    .. math::
        k[(i,j),(n,m)] := k[i,n] + k[j,m] + k[i,m] + k[j,n]

    The pairwise kernel matrix will have the same ``dtype`` as the Gram matrix
    of the given subkernel.

    Warning: may be affected by `https://github.com/numpy/numpy/issues/2396`.

    Parameters
    ----------
    pairs : ordered collection
        Pairs of indices of the elements to combine.
    subkernel : np.ndarray of shape (n, n) or Kernel instance
        Sub-kernel
    op : str, optional. (defaults to "product")
        Either "product" or "sum".
    """
    # TODO implement block-wise computation
    def __init__(self, pairs, subkernel, op="product", **kwargs):
        super(PairwiseKernel, self).__init__(pairs, **kwargs)
        self._subkernel = subkernel
        if not op in ("product", "sum"):
            raise ValueError, "invalid op '{}'".format(op)
        self._op = op

    def _compute_range_product(self, matrix, submatrix, indices1, indices2):
        for s, (i, j) in indices1:
            for t, (n, m) in indices2:
                kin, kim = submatrix[i,n], submatrix[i,m]
                kjn, kjm = submatrix[j,n], submatrix[j,m]
                matrix[s, t] = kin * kjm + kim * kjn

    def _compute_range_sum(self, matrix, submatrix, indices1, indices2):
        for s, (i, j) in indices1:
            for t, (n, m) in indices2:
                kin, kim = submatrix[i,n], submatrix[i,m]
                kjn, kjm = submatrix[j,n], submatrix[j,m]
                matrix[s, t] = kin + kjm + kim + kjn

    def _compute_all(self):
        indices = list(enumerate(self._entities))
        try:
            submatrix = self._subkernel.compute()
        except AttributeError, e:
            submatrix = self._subkernel

        matrix = np.zeros((len(self), len(self)), dtype=submatrix.dtype)
        if self._op == "product":
            self._compute_range_product(matrix, submatrix, indices, indices)
        else:
            self._compute_range_sum(matrix, submatrix, indices, indices)

        return matrix

class _TestEmpiricalKernel(object):
    _SUBMATRIX = np.array([
        [1.0, 0.5, 0.0, 0.0],
        [0.5, 1.0, 0.5, 0.0],
        [0.0, 0.5, 1.0, 0.5],
        [0.0, 0.0, 0.5, 1.0],
    ])

    def test_compute(self):
        INDICES = [1, 2]
        EXPECTED = np.array([[1.5, 1.0], [1.0, 1.5]])

        kernel = EmpiricalKernel(INDICES, self._SUBMATRIX, normalize=False)
        output = kernel.compute()
        assert (output == EXPECTED).all()

class _TestPairwiseKernel(object):
    _SUBMATRIX = np.arange(4).reshape(2, 2).astype(np.float32)
    _INDICES = list(product(range(2), range(2)))

    def _test_results(self, expected, *args, **kwargs):
        kernel = PairwiseKernel(*args, **kwargs)
        output = kernel.compute()
        assert (output == expected).all()

    def test_product(self):
        EXPECTED = np.array([
            [0,  0,  0,  2],
            [0,  2,  2,  6],
            [0,  2,  2,  6],
            [8, 12, 12, 18],
        ], dtype=np.float32)
        self._test_results(EXPECTED, self._INDICES, self._SUBMATRIX,
                           op="product", normalize=False)

    def test_sum(self):
        EXPECTED = np.array([
            [0,  2,  2,  4],
            [4,  6,  6,  8],
            [4,  6,  6,  8],
            [8, 10, 10, 12],
        ], dtype=np.float32)
        self._test_results(EXPECTED, self._INDICES, self._SUBMATRIX,
                           op="sum", normalize=False)
