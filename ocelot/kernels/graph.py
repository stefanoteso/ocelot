# -*- coding: utf8 -*-

import numpy as np
import scipy.linalg as la

from . import Kernel

class DiffusionKernel(Kernel):
    """The diffusion kernel between graph nodes [Kondor02]_.

    The diffusion kernel between two nodes of a graph :math:`G`, defined by:

    .. math::

        K(u,u') := (\\exp ( -\\beta H \\ ))_{u,u'}

    Here :math:`H` is the unnormalized Laplacian matrix of :math:`G`, and
    :math:`\\beta` is a user-provided smoothing parameter.

    Parameters
    ----------
    A : numpy.ndarray
        Adjacency matrix of the graph.
    beta : float, optional
        Non-negative smoothing parameter, defaults to 1.0.

    References
    ----------
    .. [1] Kondor and Lafferty, *Diffusion Kernels on Graphs and Other Discrete
           Input Spaces*, 2002.
    """
    def __init__(self, A, beta=1.0, **kwargs):
        if A.ndim != 2:
            raise ValueError("A is not a matrix")
        if A.shape[0] != A.shape[1]:
            raise ValueError("A is not a square matrix")
        super(DiffusionKernel, self).__init__([0] * A.shape[0], **kwargs)
        self._A, self._beta = A.astype(self.dtype), beta

    def _compute_all(self):
        A = self._A
        D = np.eye(A.shape[0], dtype=self.dtype) * np.sum(A, axis=1)
        E = la.expm(self._beta * (A - D))

        # Matrix exponentiation can give slightly asymmetric results; here we
        # fix them if they are small enough.
        asymmetry = np.linalg.norm(E.T - E, ord="fro")
        if asymmetry > 0.0:
            assert asymmetry < 1e-10, "the gods are playful today"
            E = 0.5 * (E.T + E)
        return E
