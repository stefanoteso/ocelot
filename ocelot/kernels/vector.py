# -*- coding: utf8 -*-

import numpy as np
import itertools as it

from . import Kernel

class LinearKernel(Kernel):
    """An explicit dot product kernel.

    :param entities: list of real-valued vectors.
    """
    def _compute_all(self):
        phis = np.array(self._entities)
        return np.dot(phis, phis.T)

class CorrelationKernel(Kernel):
    """An explicit dot product kernel over z-scores.

    Entities with no measurements should be set to all-zeroes.

    .. note::

        The computation of the standard deviation uses the Maximum
        Likelihood (i.e. unbiased) estimator.

    :param entities: list of real-valued vectors.
    """
    def _compute_all(self):
        phis = np.array(self._entities)
        m = np.mean(phis, axis = 1, keepdims = True)
        s = np.std(phis, axis = 1, keepdims = True)
        temp = np.seterr(divide = 'ignore', invalid = 'ignore')
        z = (phis - m) / s
        np.seterr(**temp)
        z[np.isnan(z)] = 0.0
        return np.dot(z, z.T)

class SparseLinearKernel(Kernel):
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

class SetKernel(Kernel):
    """A generic set kernel.

    :param entities: list of sets.
    """
    def _compute_all(self):
        num = len(self)
        matrix = np.zeros((num, num))
        for i in xrange(num):
            entity_i = self._entities[i]
            for j in xrange(i + 1):
                entity_j = self._entities[j]
                matrix[i,j] = matrix[j,i] = float(len(entity_i & entity_j))
        return matrix

class ColocalizationKernel(Kernel):
    """An exponential kernel over genetic contexts.

    :param contexts: a sequence of tuples of the form ``(chromosome, position)``.
    """
    def __init__(self, contexts, *args, **kwargs):
        self._gamma = kwargs.get("gamma", 1.0)
        super(ColocalizationKernel, self).__init__(contexts, *args, **kwargs)

    def _compute_all(self):
        # Group the contexts by chromosome
        chromosome_to_contexts = {}
        for i, (chromosome, pos) in enumerate(self._entities):
            if not chromosome in chromosome_to_contexts:
                chromosome_to_contexts[chromosome] = []
            chromosome_to_contexts[chromosome].append((i, pos))

        matrix = np.zeros((len(self), len(self)))
        for contexts in chromosome_to_contexts.itervalues():
            # Figure out the maximum distance between genes
            max_d = None
            for (i, pos_i), (j, pos_j) in it.product(contexts, contexts):
                if i <= j:
                    continue
                d = np.abs(pos_i - pos_j)
                if d > max_d or max_d == None:
                    max_d = d
            # Compute the kernel matrix
            for (i, pos_i), (j, pos_j) in it.product(contexts, contexts):
                if i < j:
                    continue
                matrix[i,j] = \
                matrix[j,i] = \
                     np.exp(-self._gamma * (np.abs(pos_i - pos_j) / max_d))
        return matrix

class _TestLinearKernel(object):
    def test_result(self):
        INPUTS = (
            ((np.array([0, 0]), np.array([0, 0])), np.array([[0, 0], [0, 0]])),
            ((np.array([1, 0]), np.array([0, 1])), np.array([[1, 0], [0, 1]])),
            ((np.array([1, 0]), np.array([1, 0])), np.array([[1, 1], [1, 1]])),
        )
        for phis, expected in INPUTS:
            kernel = LinearKernel(phis, do_normalize = False)
            output = kernel.compute()
            assert (output == expected).all()

class _TestCorrelationKernel(object):
    def test_result(self):
        INPUTS = (
            ((np.array([0, 0]), np.array([0, 0])), np.array([[0, 0], [0, 0]])),
            ((np.array([1, 0]), np.array([0, 1])), np.array([[2, -2], [-2, 2]])),
            ((np.array([1, 0]), np.array([1, 0])), np.array([[2, 2], [2, 2]])),
        )
        for phis, expected in INPUTS:
            kernel = CorrelationKernel(phis, do_normalize = False)
            output = kernel.compute()
            assert (output == expected).all()

class _TestSetKernel(object):
    def test_results(self):
        INPUTS = (
            ((set(), set()), np.array([[0, 0], [0, 0]])),
            ((set([0]), set([1])), np.array([[1, 0], [0, 1]])),
            ((set([0]), set([0])), np.array([[1, 1], [1, 1]])),
            ((set(["a", "b"]), set(["a", "b"])), np.array([[2, 2], [2, 2]])),
        )
        for phis, expected in INPUTS:
            kernel = SetKernel(phis, do_normalize = False)
            output = kernel.compute()
            assert (output == expected).all()
