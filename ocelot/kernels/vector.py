# -*- coding: utf8 -*-

import numpy as np

from .base import _Kernel

class LinearKernel(_Kernel):
    """An explicit dot product kernel.

    :param entities: list of real-valued vectors.
    """
    def _compute_all(self):
        num = len(self)
        matrix = np.zeros((num, num))
        for i in xrange(num):
            phi_i = self._entities[i]
            for j in xrange(i + 1):
                phi_j = self._entities[j]
                matrix[i,j] = matrix[j,i] = np.dot(phi_i.T, phi_j)
        return matrix

class CorrelationKernel(_Kernel):
    """A linear kernel over z-scores (computed autmatically).

    :param entities: list of real-valued vectors.
    """
    def _compute_all(self):
        for i, vector in enumerate(self._entities):
            m = np.mean(vector)
            s = np.std(vector)
            if s == 0.0:
                assert m == 0.0
                # missing expression levels for this entry
                s = 1.0
            self._entities[i] = (vector - m) / s
        num = len(self)
        matrix = np.zeros((num, num))
        for i in xrange(num):
            phi_i = self._entities[i]
            if not any(phi_i):
                # missing expression levels for this entry
                matrix[i,i] = 1.0
                continue
            for j in xrange(i + 1):
                phi_j = self._entities[j]
                matrix[i,j] = matrix[j,i] = np.dot(phi_i.T, phi_j)
        return matrix

class SparseLinearKernel(_Kernel):
    """A sparse linear kernel.

    :param entities: list of dictionaries whose keys are the indices and the
        corresponding values are the weights associated to the indices.
    """
    def _compute_all(self):
        num = len(self)
        matrix = np.zeros((num, num), dtype=np.float64)
        for i in xrange(num):
            entity_i = self._entities[i]
            for j in xrange(i + 1):
                entity_j = self._entities[j]
                if i == j and len(entity_i) == 0 and len(entity_j) == 0:
                    dp = 1.0
                else:
                    common_keys = set(entity_i.keys()) & set(entity_j.keys())
                    dp = sum([entity_i[k]*entity_j[k] for k in common_keys])
                matrix[i,j] = matrix[j,i] = dp
        return matrix

class SetKernel(_Kernel):
    """A generic set kernel.

    :param entities: list of sets.
    """
    def _compute_all(self):
        num = len(self)
        matrix = np.zeros((num, num), dtype=np.float64)
        for i in xrange(num):
            entity_i = self._entities[i]
            for j in xrange(i + 1):
                entity_j = self._entities[j]
                if i == j and len(entity_i) == 0 and len(entity_j) == 0:
                    dp = 1.0
                else:
                    dp = float(len(entity_i & entity_j))
                matrix[i,j] = matrix[j,i] = dp
        return matrix
