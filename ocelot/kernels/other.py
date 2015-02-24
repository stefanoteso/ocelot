# -*- coding: utf8 -*-

from .base import _Kernel
import numpy as np

class EmpiricalKernelMap(_Kernel):
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

class PairwiseKernel(_Kernel):
    """The pairwise kernel.

    It computes a kernel between pairs of elements as follows:

    .. math::
        k[(i,j),(n,m)] := k[i,n]k[j,m] + k[i,m]k[j,n]

    or:

    .. math::
        k[(i,j),(n,m)] := k[i,n] + k[j,m] + k[i,m] + k[j,n]

    :param indices: indices of the elements to combine.
    :param subkernel: a `Kernel` instance.
    :param op: either ``"tensor"`` or ``"sum"``.
    """
    def __init__(self, pairs, subkernel, op = "tensor", **kwargs):
        super(PairwiseKernel, self).__init__(pairs, **kwargs)
        self.subkernel = subkernel
        if op == "tensor":
            self.op = lambda kin, kjm, kim, kjn: kin * kjm + kim * kjn
        elif op == "sum":
            self.op = lambda kin, kjm, kim, kjn: kin + kjm + kim + kjn
        else:
            raise ValueError, "invalid op '{}'".format(op)
    def _compute_all(self):
        submatrix = self.subkernel.compute()
        matrix = np.zeros((len(self), len(self)))
        for out_i, (i, j) in enumerate(self._entities):
            for out_j, (n, m) in enumerate(self._entities):
                kin, kim = submatrix[i,n], submatrix[i,m]
                kjn, kjm = submatrix[j,n], submatrix[j,m]
                matrix[out_i, out_j] = self.op(kin, kjm, kim, kjn)
        return matrix
