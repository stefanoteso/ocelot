# -*- coding: utf8 -*-

import numpy as np
from ocelot.services import AMINOACIDS, BACKGROUND_AA_FREQ, PSSM

from .base import _Kernel

class _RecursivePrefixStringKernel(_Kernel):
    """Base prefix-tree string kernel implemented using recursion.

    Used to implement the spectrum, mismatch and profile kernels. Adapted from
    `<http://cbio.mskcc.org/leslielab/software/string_kernels.html>`_.

    Note that more efficient, non-recursive implementations are possible. See
    for instance the `fastprof <https://rostlab.org/owiki/index.php/Fastprofkernel>`_
    kernel.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_instances_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)

    .. todo::

        The current formulation completely throws away information about
        *l*-mers, where :math:`l < k`. We should avoid this, and compute all
        kernel matrices for :math:`l \leq k` in a single pass. In the original
        code this is done by aggregating instances at each node, rather than
        filtering them out like we do here, and storing the number of hops
        survived in the instances themselves, then adding a check in _update.
        However this is **very** slow.
    """
    def __init__(self, strings, **kwargs):
        super(_RecursivePrefixStringKernel, self).__init__(strings, **kwargs)
        self._k = kwargs.get("k", 1)
        if not self._k >= 1:
            raise ValueError("k must be >= 1")
        self._alphabet = kwargs.get("alphabet", AMINOACIDS)
        self._min_instances_per_node = kwargs.get("min_instances_per_node", 1)
        if not self._min_instances_per_node >= 1:
            raise ValueError("min_instances_per_node must be >= 1")

    def _update(self, matrix, instances):
        """Called when a leaf of the prefix tree is visited. Updates the kernel
        matrix based on the instances that managed to not-be filtered out."""
        raise NotImplementedError

    def _recurse(self, matrix, instances, depth, prefix):
        """Depth-first traversal of the prefix-tree."""
        if depth == self._k:
            self._update(matrix, instances)
        else:
            for symbol in self._alphabet:
                # Filter out all instances that do not pass the check at the
                # current node 
                filtered_instances = []
                for instance in instances:
                    new_instance, check = self._check_instance(instance, symbol, depth)
                    if check:
                        filtered_instances.append(new_instance)
                # If there are no instances left, there is nothing else to do
                if len(filtered_instances) < self._min_instances_per_node:
                    continue
                self._recurse(matrix, filtered_instances, depth + 1, prefix + symbol)
                del filtered_instances

    def _get_instances(self):
        """Turns the input strings into instances, i.e. k-mers plus additional
        information."""
        instances = []
        for i, string in enumerate(self._entities):
            if len(string) < self._k:
                raise ValueError("string shorter than k")
            instances.extend([self._substring_to_instance(i, offset)
                              for offset in xrange(len(string) - self._k + 1)])
        return instances

    def _compute_all(self):
        matrix = np.zeros((len(self), len(self)))
        self._instances = self._get_instances()
        self._recurse(matrix, self._instances, 0, "")
        return matrix

class SpectrumKernel(_RecursivePrefixStringKernel):
    """The spectrum kernel for strings [Leslie02a]_.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_instances_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def _update(self, matrix, instances):
        for i, _ in instances:
            for j, _ in instances:
                matrix[i,j] += 1

    def _check_instance(self, instance, symbol, depth):
        i, offset = instance
        string_symbol = self._entities[i][offset + depth]
        return instance, string_symbol == symbol

    def _substring_to_instance(self, i, offset):
        return (i, offset)

class MismatchKernel(_RecursivePrefixStringKernel):
    """The mismatch kernel for strings [Leslie02b]_.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param m: maximum number of mismatches (default: 0).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_instances_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def __init__(self, strings, **kwargs):
        super(MismatchKernel, self).__init__(strings, **kwargs)
        self._m = kwargs.get("m", 1)
        if not self._m >= 0:
            raise ValueError("m must be >= 0")

    def _update(self, matrix, instances):
        for i, _, _ in instances:
            for j, _, _ in instances:
                matrix[i,j] += 1

    def _check_instance(self, instance, symbol, depth):
        i, offset, num_mismatches = instance
        string_symbol = self._entities[i][offset + depth]
        if string_symbol != symbol:
            num_mismatches += 1
        new_instance = (i, offset, num_mismatches)
        return new_instance, num_mismatches <= self._m

    def _substring_to_instance(self, i, offset):
        return (i, offset, 0)

class ProfileKernel(_RecursivePrefixStringKernel):
    """Implementation of the PSSM-based string kernel [Kuang04]_.

    :param pssms: list of PSSM matrices.
    :param k: kmer size, inclusive (default: 1).
    :param threshold: threshold mutation probability to count as a hit (default: 6.0)
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param gamma: amount of smoothing-by-prior (default: 0.8).
    :param prior: list of priors on amino acid probabilities (default: BACKGROUND_AA_FREQ)
    :param min_instances_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def __init__(self, pssms, **kwargs):
        super(ProfileKernel, self).__init__(pssms, **kwargs)
        self._threshold = kwargs.get("threshold", 6.0)
        self._gamma = kwargs.get("gamma", 0.8)
        if not (0.0 <= self._gamma <= 1.0):
            raise ValueError("gamma must be in [0.0, 1.0]")
        self._prior = kwargs.get("prior", BACKGROUND_AA_FREQ)

        # We convert the PSSM transition probabilities into transition scores
        # by (i) smoothing the observed transition probabilities given by the
        # PSSM by some background frequency prior, and (ii) taking the negative
        # log of the resulting mixture.
        self._scores = map(self._to_scores, pssms)

    def _to_scores(self, pssm):
        scores = []
        for i, (in_aminoacid, transition_prob) in enumerate(pssm):
            row = []
            for j, out_aminoacid in enumerate(self._alphabet):
                p1 = transition_prob[j]
                p2 = self._prior[out_aminoacid]
                p = self._gamma * p1 + (1 - self._gamma) * p2
                assert 0.0 <= p <= 1.0
                if p == 0.0:
                    # XXX some huge value
                    row.append(1e13)
                else:
                    row.append(-np.log(p))
                assert row[-1] >= 0.0
            scores.append(row)
        return scores

    def _update(self, matrix, instances):
        for i, _, _ in instances:
            for j, _, _ in instances:
                matrix[i, j] += 1

    def _check_instance(self, instance, symbol, depth):
        i, offset, score = instance
        j = self._alphabet.index(symbol)
        transition_score = score + self._scores[i][offset + depth][j]
        new_instance = (i, offset, transition_score)
        # Recall that more likely means (i) probability closer to 1, but (ii)
        # score (= neg log of probability) closer to 0. So to capture the most
        # likely mutant k-mers we want the total score not to be increase too
        # much.
        return new_instance, transition_score <= self._threshold

    def _substring_to_instance(self, i, offset):
        return (i, offset, 0.0)
