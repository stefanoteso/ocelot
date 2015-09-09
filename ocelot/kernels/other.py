# -*- coding: utf8 -*-

import numpy as np
import itertools as it

from . import Kernel, DummyKernel

class EmpiricalKernelMap(Kernel):
    """Empirical kernel map over an underlying kernel[1].

    *References*

    [1] Scholkopf et al., "Input Space Versus Feature Space in Kernel-Based
        Methods", 1999.
    """
    def __init__(self, *args, **kwargs):
        import ocelot.features
        super(EmpiricalKernelMap, self).__init__(*args, **kwargs)
        self.features = features.EmpiricalFeatures(args[0], range(len(self)))
    def _compute_all(self):
        LinearKernel(self.features, range(len(self)))

class PairwiseKernel(Kernel):
    """The pairwise kernel.

    It computes a kernel between pairs of elements as follows:

    .. math::
        k[(i,j),(n,m)] := k[i,n]k[j,m] + k[i,m]k[j,n]

    or:

    .. math::
        k[(i,j),(n,m)] := k[i,n] + k[j,m] + k[i,m] + k[j,n]

    The pairwise kernel matrix will have the same ``dtype`` as the Gram matrix
    of the given subkernel.

    .. warning::

        Make sure not to be affected by `https://github.com/numpy/numpy/issues/2396`.

    .. todo::

        Implement block-by-block computation.

    :param indices: indices of the elements to combine.
    :param subkernel: a `Kernel` instance or an ``numpy.ndarray``.
    :param op: either ``"product"`` or ``"sum"`` (default: ``"product"``).
    """
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
        assert len(indices) == len(set(indices)), "repeated indices found"
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

class _TestPairwiseKernel(object):
    _SUBMATRIX = np.arange(4).reshape(2, 2).astype(np.float32)
    _INDICES = list(it.product(range(2), range(2)))

    def _test_results(self, expected, *args, **kwargs):
        kwargs["do_normalize"] = False
        kernel = PairwiseKernel(*args, **kwargs)
        output = kernel.compute()
        assert (output == expected).all()

    def test_tensor_product(self):
        EXPECTED = np.array([
            [0,  0,  0,  2],
            [0,  2,  2,  6],
            [0,  2,  2,  6],
            [8, 12, 12, 18],
        ], dtype=np.float32)
        self._test_results(EXPECTED, self._INDICES, self._SUBMATRIX, op="product")

    def test_tensor_sum(self):
        EXPECTED = np.array([
            [0,  2,  2,  4],
            [4,  6,  6,  8],
            [4,  6,  6,  8],
            [8, 10, 10, 12],
        ], dtype=np.float32)
        self._test_results(EXPECTED, self._INDICES, self._SUBMATRIX, op="sum")
