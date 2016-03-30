# -*- coding: utf8 -*-

import numpy as np
import cPickle as pickle
from ocelot.utils import ispsd

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
except ImportError:
    pass

_EPS = 1e-10

def kernnorm(matrix):
    """Computes the normalized kernel matrix.

    NOTE: zero-valued diagonal elements are gracefully dealt with; but
    negative diagonal elements will lead to NaNs (and asymmetric matrices
    due to the NaN comparison rules).

    .. math::

        \\hat{k}_{ij} = k_{ij} / \\sqrt{k_ii k_jj}

    Parameters
    ----------
    matrix : np.ndarray of shape (n, n)
        The Gram matrix.
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
    if matrix2_or_y.ndim == 1:
        matrix2 = np.outer(matrix2_or_y, matrix2_or_y).ravel()
    else:
        matrix2 = matrix2_or_y.ravel()

    if indices is not None:
        matrix1 = matrix1[np.ix_(indices, indices)]
        matrix2 = matrix2[np.ix_(indices, indices)]

    matrix1, matrix2 = matrix1.ravel(), matrix2.ravel()
    dot11 = np.dot(matrix1, matrix1)
    dot12 = np.dot(matrix1, matrix2)
    dot22 = np.dot(matrix2, matrix2)
    return dot12 / (np.sqrt(dot11 * dot22) + _EPS)

def draw_kernel(matrix, path):
    """Draws a Gram matrix."""
    figsize=(matrix.shape[0] // 96, matrix.shape[1] // 96)
    fig = plt.figure(figsize=figsize, dpi=96)
    ax = fig.add_subplot(111)
    ax.matshow(matrix, interpolation="nearest", cmap=cm.OrRd)
    fig.savefig(path, bbox_inches="tight", pad_inches=0)

class Kernel(object):
    """Base kernel class.

    Parameters
    ----------
    entities : collection
        Ordered collection of arbitrary elements.
    normalize : bool, optional. (defaults to False)
        Whether to normalize the kernel.
    fixup : False or pair of floats. (defaults to False)
        Whether to fixup the kernel, and by how much.
    dtype : np.dtype, optional. (defaults to None)
        WRITEME
    """
    def __init__(self, entities, normalize=False, fixup=False, dtype=None):
        if len(entities) < 2:
            raise ValueError("kernel need at least two entities.")
        self._entities = entities
        self._normalize = normalize
        self._fixup = fixup
        self.dtype = dtype

    def ispsd(self, tol=1e-6):
        return ispsd(self.compute(), tol=tol)

    def normalize(self):
        self._matrix = kernnorm(self.compute())

    def fixup(self, asymm_tol=1e-6, nonpsd_tol=1e-6):
        matrix = self.compute()

        matrix = 0.5 * (matrix + matrix.T)
        ls = np.linalg.eigvalsh(matrix)
        if ls[0] < 0.0:
            if ls[0] >= nonpsd_tol:
                raise ValueError("matrix is too non-PSD '{}'".format(ls[0]))
            matrix += np.eye(matrix.shape[0]) * -10.0 * ls[0]

        self._matrix = matrix

    def alignment(self, matrix_or_y, indices=None):
        return kernalign(self.compute(), matrix_or_y, indices=indices)

    def draw(self, path):
        draw_kernel(self.compute(), path)

    def __len__(self):
        """Returns the number of entities."""
        return len(self._entities)

    def compute(self):
        """Computes the Gram matrix.

        Calls `self._compute_all`.
        """
        if not hasattr(self, "_matrix"):
            self._matrix = np.array(self._compute_all())
            if self._normalize:
                self.normalize()
            if self._fixup != False:
                self.fixup(asymm_tol=self._fixup[0],
                           nonpsd_tol=self._fixup[1])
        return self._matrix

class DummyKernel(Kernel):
    """A wrapper around ``np.ndarray``'s and files.

    Parameters
    ----------
    matrix_or_path : numpy.ndarray or str
        Either a Gram matrix or a path to a file.
    All other arguments are passed to the ``Kernel`` constructor.
    """
    def __init__(self, matrix_or_path, **kwargs):
        matrix = self._load(matrix_or_path)
        if matrix is None:
            raise ValueError("invalid matrix_or_path")

        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("matrix is not square '{}'".format(matrix.shape))

        super(DummyKernel, self).__init__(range(matrix.shape[0]),
                                          **kwargs)

        self._matrix = matrix.astype(self.dtype)
        if self._normalize:
            self.normalize()
        if self._fixup != False:
            self.fixup(asymm_tol=self._fixup[0],
                       nonpsd_tol=self._fixup[1])

    @staticmethod
    def _load(matrix_or_path):
        matrix = None
        try:
            matrix_or_path.shape
            matrix = matrix_or_path
        except:
            pass
        try:
            if matrix is None:
                matrix = np.loadtxt(matrix_or_path)
        except:
            pass
        try:
            if matrix is None:
                with open(matrix_or_path, "rb") as fp:
                    matrix = pickle.load(fp)
        except:
            pass
        return matrix

    def _compute_all(self):
        return self._matrix

class _TestKernel(object):
    def _test(self, matrix, expected, normalize=True):
        kernel = DummyKernel(np.array(matrix), normalize=normalize)
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
