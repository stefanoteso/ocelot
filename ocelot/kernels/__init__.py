# -*- coding: utf8 -*-

import numpy as np

_EPS = 1e-10

def kernalign(matrix1, matrix2_or_y, indices=None):
    """Computes the kernel-kernel or kernel-target alignment [1]_:

    .. math::
        <matrix1, matrix2> / \\sqrt{ <matrix1, matrix1> <matrix2, matrix2> }

    Parameters
    ----------
    matrix1 : np.ndarray of shape (n, n)
        The Gram matrix.
    matrix2_or_target : np.ndarray of shape (n, n) or (n,)
        Either a second Gram matrix or a target vector. If the shape is (n,),
        then matrix2 is the outer product :math:`target target^\\top`.
    indices : ordered collection
        Indices that define the empirical sample.

    Returns
    -------
    alignment : float
        The alignment.

    References
    ----------
    .. [1] N. Cristianini et al. *On Kernel-Target Alignment*, 2001.
    """
    vec1 = matrix1.ravel()
    if matrix2_or_y.ndim == 1:
        vec2 = np.outer(matrix2_or_y, matrix2_or_y).ravel()
    else:
        vec2 = matrix2_or_y.ravel()
    dot11 = np.dot(vec1, vec1)
    dot12 = np.dot(vec1, vec2)
    dot22 = np.dot(vec2, vec2)
    return dot12 / (np.sqrt(dot11 * dot22) + _EPS)

class Kernel(object):
    """Base kernel class.

    Parameters
    ----------
    entities : collection
        Ordered collection of arbitrary elements.
    do_normalize : boo, optional (defaults to True)
        Whether to normalize the Gram matrix.
    """
    def __init__(self, entities, *args, **kwargs):
        self._entities = entities
        if not (len(self._entities) >= 2):
            raise ValueError("kernels need at least two entities.")
        self._do_normalize = kwargs.get("do_normalize", True)
        self.dtype = kwargs.get("dtype", np.float32)
        self._matrix = None

    def __len__(self):
        """Returns the number of entities."""
        return len(self._entities)

    @staticmethod
    def _normalize(matrix):
        """Returns the normalized kernel matrix.

        NOTE: zero-valued diagonal elements are gracefully dealt with; but
        negative diagonal elements will lead to NaNs (and asymmetric matrices
        due to the NaN comparison rules).

        .. math::

            `\hat{k}_{ij} = k_{ij} / \sqrt{k_ii k_jj}`
        """
        diag = matrix.diagonal()
        num_zero_diag_entries = np.where(diag == 0)[0].shape[0]
        if num_zero_diag_entries:
            print "Warning: found zero-valued diagonal entries during normalization, will fix up INF to zeros."
        z = 1.0 / np.sqrt(diag)
        # We may get infinities if the diagonal is zero; in that case we
        # *assume* that the corresponding non-diagonal entries are also zero,
        # which implies that it is safe to turn the infinities to 1, as doing
        # so retains the zeros on the non-diagonal entries.
        z[np.isinf(z)] = 1.0
        for i, row in enumerate(matrix):
            matrix[i] = row * z * z[i]
        return matrix

    def compute(self):
        """Computes the Gram matrix.

        Relies on `self._compute_all`.
        """
        if self._matrix is None:
            self._matrix = np.array(self._compute_all())
            if self._do_normalize:
                self._matrix = self._normalize(self._matrix)
        return self._matrix

    def is_psd(self):
        """Checks whether the Gram matrix is positive semi-definite."""
        matrix = self.compute()
        is_symmetric = (np.abs(matrix.T - matrix) < 1e-6).all()
        is_semi_psd = np.all(np.linalg.eigvalsh(matrix) >= 0)
        return is_symmetric and is_semi_psd

    def check_and_fixup(self, threshold):
        """Checks whether the Gram matrix is positive semi-definite and
        preconditions it if it is slightly non-PSD."""
        matrix = self.compute()
        assert (np.abs(matrix.T - matrix) < 1e-6).all(), "not symmetric!"
        matrix = 0.5 * (matrix + matrix.T)

        ls = np.linalg.eigvalsh(matrix)
        if ls[0] < 0.0:
            assert ls[0] < threshold, "matrix is too non-PSD: minimum eigenvalue is '{}'".format(ls[0])
            print "preconditioning by 10*(ls[0] = '{}')".format(ls[0])
            matrix += np.identity(matrix.shape[0]) * -10.0 * ls[0]
            eigvals = np.linalg.eigvalsh(matrix)
            assert eigvals[0] > -threshold, "eigenvalues are too non-negative, even after preconditioning: {}".format(eigvals)

    def draw(self, path):
        """Draws the Gram matrix to a file."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm
        except ImportError:
            print "matplotlib is required; can not draw"
            return
        try:
            fig = plt.figure(figsize = (len(self) / 96, len(self) / 96), dpi = 96)
            ax = fig.add_subplot(1, 1, 1)
            matrix = self.compute()
            ax.matshow(matrix, interpolation = "nearest", cmap = cm.OrRd)
            fig.savefig(path, bbox_inches="tight", pad_inches=0)
        except Exception, e:
            print "failed to draw kernel; skipping"

class DummyKernel(Kernel):
    """A wrapper around ``np.ndarray``'s and files.

    Parameters
    ----------
    arg : numpy.ndarray or str
        Either a Gram matrix or a path to a file.
    num : int
        Number of elements in the Gram matrix.
    check_psd : bool
        Whether to raise an exception if the wrapped kernel is not PSD.
    """
    def __init__(self, arg, **kwargs):
        arg = arg
        num = kwargs.get("num")
        check_psd = kwargs.get("check_psd")
        if isinstance(arg, str):
            matrix = np.loadtxt(arg)
        else:
            matrix = arg
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("matrix is not symmetric '{}'".format(matrix.shape))
        if not num is None and matrix.shape != (num, num):
            raise ValueError("matrix has invalid shape '{}', was expecting ({},{})".format(matrix.shape, num, num))
        super(DummyKernel, self).__init__(range(matrix.shape[0]),
                                          **kwargs)
        self._matrix = matrix
        if check_psd and not self.is_psd():
            raise ValueError("matrix is not PSD")
        if self._do_normalize:
            self._matrix = self._normalize(self._matrix)

class _TestKernel(object):
    def _test(self, matrix, expected, do_normalize=True):
        kernel = DummyKernel(np.array(matrix), do_normalize=do_normalize)
        output = kernel.compute()
        assert (np.abs(output - expected) < 1e-6).all()

    def test_normalization(self):
        MATRIX = np.array([
            [2, 1, 0],
            [1, 2, 1],
            [0, 1, 2],
        ], dtype=np.float64)
        EXPECTED = np.array([
            [1,   0.5,   0],
            [0.5,   1, 0.5],
            [0,   0.5,   1],
        ], dtype=np.float64)
        self._test(MATRIX, EXPECTED)

    def test_normalization_zero_diag(self):
        # zero diagonal elements are treated as ones
        S = 1.0 / np.sqrt(2)
        MATRIX = np.array([
            [2, 1, 0],
            [1, 0, 1],
            [0, 1, 2],
        ], dtype=np.float64)
        EXPECTED = np.array([
            [1, S, 0],
            [S, 0, S],
            [0, S, 1],
        ], dtype=np.float64)
        self._test(MATRIX, EXPECTED)

from .vector import *
from .string import *
from .graph import *
from .other import *
