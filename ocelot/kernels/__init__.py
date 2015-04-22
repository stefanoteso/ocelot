# -*- coding: utf8 -*-

import numpy as np

from ..services import _cls, AMINOACIDS, BACKGROUND_AA_FREQ

class Kernel(object):
    """Base kernel class.

    :param entities: list of arbitrary data entries.
    :param do_normalize: normalize the Gram matrix.
    """
    def __init__(self, entities, **kwargs):
        self._entities = entities
        if not (len(self._entities) >= 2):
            raise ValueError("kernels need at least two entities.")
        self._do_normalize = kwargs.get("do_normalize", True)
        self._matrix = None

    def __len__(self):
        return len(self._entities)

    @staticmethod
    def _normalize(matrix):
        """Returns the normalized kernel matrix.

        .. math::

            `\hat{k}_{ij} = k_{ij} / \sqrt{k_ii k_jj}`
        """
        ones = np.ones((1, matrix.shape[1]))
        invd = np.divide(ones, np.sqrt(matrix.diagonal()))
        return np.multiply(matrix, np.dot(invd.T, invd))

    def compute(self):
        """Computes the kernel matrix.

        Relies on `self._compute_all`.
        """
        if self._matrix == None:
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
            assert np.linalg.eigvalsh(matrix)[0] >= 0, "the gods are playful today"

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
        if num != None and matrix.shape != (num, num):
            raise ValueError("matrix has invalid shape '{}', was expecting ({},{})".format(matrix.shape, num, num))
        super(DummyKernel, self).__init__(range(matrix.shape[0]),
                                          **kwargs)
        self._matrix = matrix
        if check_psd and not self.is_psd():
            raise ValueError("matrix is not PSD")
        if self._do_normalize:
            self._matrix = self._normalize(self._matrix)

from .vector import *
from .string import *
from .graph import *
from .other import *
