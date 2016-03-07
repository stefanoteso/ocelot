# -*- coding: utf8 -*-

import numpy as np

from ..services import _cls

class Kernel(object):
    """Base kernel class.

    :param entities: list of arbitrary data entries.
    :param do_normalize: normalize the Gram matrix.
    """
    def __init__(self, entities, *args, **kwargs):
        self._entities = entities
        if not (len(self._entities) >= 2):
            raise ValueError("kernels need at least two entities.")
        self._do_normalize = kwargs.get("do_normalize", True)
        self.dtype = kwargs.get("dtype", np.float32)
        self._matrix = None

    def __len__(self):
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
        """Computes the kernel matrix.

        Relies on `self._compute_all`.
        """
        if self._matrix is None:
            self._matrix = self._compute_all()
            if self._do_normalize:
                self._matrix = self._normalize(self._matrix)
        return self._matrix

    def is_psd(self):
        matrix = self.compute()
        is_symmetric = (matrix.T == matrix).all()
        is_semi_psd = np.all(np.linalg.eigvalsh(matrix) >= 0)
        return is_symmetric and is_semi_psd

    def check_and_fixup(self, threshold):
        matrix = self.compute()
        assert (matrix.T == matrix).all(), "not symmetric!"

        ls = np.linalg.eigvalsh(matrix)
        if ls[0] < 0.0:
            assert ls[0] < threshold, "matrix is too non-PSD: minimum eigenvalue is '{}'".format(ls[0])
            print _cls(self), ": preconditioning by 10*(ls[0] = '{}')".format(ls[0])
            matrix += np.identity(matrix.shape[0]) * -10.0 * ls[0]
            eigvals = np.linalg.eigvalsh(matrix)
            assert eigvals[0] > -threshold, "eigenvalues are too non-negative, even after preconditioning: {}".format(eigvals)

    def draw(self, path):
        """Draws the kernel matrix to a file."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm
        except:
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

    def save(self, path):
        np.savetxt(path, self.compute())

class DummyKernel(Kernel):
    """A wrapper around ``np.ndarray``'s and files.

    :param arg: either a ``numpy.ndarray`` or a path to a file.
    :param num: number of elements in the wrapper kernel.
    :param check_psd: whether to raise an exception if the wrapped kernel is not PSD.
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
    def _test(matrix, expected):
        kernel = DummyKernel(np.array(matrix), do_normalize=True)
        output = kernel.compute()
        assert (output == expected).all()

    def test_normalization(self):
        MATRIX = np.array([
            [2, 1, 0],
            [1, 2, 1],
            [0, 1, 2],
        ])
        EXPECTED = np.array([
            [1, 1/2., 0],
            [1/2., 1, 1/2.],
            [0, 1/2., 1],
        ])
        self._test(MATRIX, EXPECTED)

    def test_normalization_weird(self):
        MATRIX = np.array([
            [2, 1, 0],
            [1, 0, 1],
            [0, 2, 2],
        ])
        EXPECTED = np.array([
            [1, 1, 0],
            [1, 1, 1],
            [0, 1, 1],
        ])
        self._test(MATRIX, EXPECTED)

from .vector import *
from .string import *
from .graph import *
from .other import *
