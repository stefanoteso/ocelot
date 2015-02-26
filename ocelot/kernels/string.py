# -*- coding: utf8 -*-

import numpy as np
from ocelot.services import AMINOACIDS

from .base import _Kernel

class _RecursivePrefixStringKernel(_Kernel):
    """Base prefix-tree string kernel implemented using recursion.

    Used to implement the spectrum, mismatch and profile kernels. Adapted from
    `<http://cbio.mskcc.org/leslielab/software/string_kernels.html>`_.

    The matrix update mechanism has been taken from
    `fastprof <https://rostlab.org/owiki/index.php/Fastprofkernel>`_.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_survivors_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)

    .. todo::

        The current formulation completely throws away information about
        *l*-mers, where :math:`l < k`. We should avoid this, and compute all
        kernel matrices for :math:`l \leq k` in a single pass. In the original
        code this is done by aggregating instances at each node, rather than
        filtering them out like we do here, and storing the number of hops
        survived in the instances themselves, then adding a check in _update.
        However this is **very** slow.

        An alternative is to update the matrix when all child nodes have
        been considered.
    """
    def __init__(self, strings, **kwargs):
        super(_RecursivePrefixStringKernel, self).__init__(strings, **kwargs)
        self._k = kwargs.get("k", 1)
        if not self._k >= 1:
            raise ValueError("k must be >= 1")
        self._alphabet = kwargs.get("alphabet", AMINOACIDS)
        self._min_survivors_per_node = kwargs.get("min_survivors_per_node", 1)
        if not self._min_survivors_per_node >= 1:
            raise ValueError("min_survivors_per_node must be >= 1")
        self._survivors = []
        # Precompute a couple things
        self._index_symbol_pairs = list(enumerate(self._alphabet))

    def _add_survivors(self, instances):
        for instance_i in instances:
            i = instance_i[0]
            for instance_j in instances:
                j = instance_j[0]
                self._matrix[i,j] += 1

#    def _add_survivors(self, instances, leaf):
#        """Called upon visiting a leaf of the prefix tree. Appends the list of
#        survivors with the surviving instances."""
#        raise NotImplementedError
#
#        for i in xrange(self._entities):
#            self._survivors.append((i, leaf, len(filter(lambda instance: instance[0] == i, instances))))
#        if len(self._survivors) >= self._max_survivors_per_flush:
#            self._survivors = []

    def _recurse(self, instances, depth, prefix):
        """Depth-first traversal of the prefix-tree."""
        if depth == self._k:
            self._add_survivors(instances)
        else:
            for s, symbol in self._index_symbol_pairs:
                # Filter out all instances that do not pass the check at the
                # current node 
                survivors = []
                for instance in instances:
                    new_instance, check = self._check_instance(instance, s, symbol, depth)
                    if check:
                        survivors.append(new_instance)
                # If there are no instances left, there is nothing else to do
                if len(survivors) < self._min_survivors_per_node:
                    continue
                # Process the child nodes
                self._recurse(survivors, depth + 1, prefix + symbol)
                # Get rid of the list of survivors
                del survivors

    def _to_instance(self, i, offset):
        """Given an element index and an offest, returns an instance (i.e.
        the k-mer indices optionally associated to additional information).

        **Note**: the first element of an instance **must** be the index of
        the pssm it refers to! The ``_add_survivors`` method relies on it.
        """
        raise NotImplementedError

    def _get_instances(self):
        """Turns the input strings into instances, i.e. k-mers plus additional
        information.
        """
        instances = []
        for i, string in enumerate(self._entities):
            if len(string) < self._k:
                raise ValueError("string shorter than k")
            instances.extend([self._to_instance(i, offset)
                              for offset in xrange(len(string) - self._k + 1)])
        return instances

    def _compute_all(self):
        self._matrix = np.zeros((len(self), len(self)))
        self._instances = self._get_instances()
        self._recurse(self._instances, 0, "")
        self._add_survivors([], -1)
        return self._matrix

class SpectrumKernel(_RecursivePrefixStringKernel):
    """The spectrum kernel for strings [Leslie02a]_.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_survivors_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def _check_instance(self, instance, _, symbol, depth):
        i, offset = instance
        string_symbol = self._entities[i][offset + depth]
        return instance, string_symbol == symbol

    def _to_instance(self, i, offset):
        return (i, offset)

class MismatchKernel(_RecursivePrefixStringKernel):
    """The mismatch kernel for strings [Leslie02b]_.

    :param strings: a list of strings.
    :param k: kmer size, inclusive (default: 1).
    :param m: maximum number of mismatches (default: 0).
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_survivors_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def __init__(self, strings, **kwargs):
        super(MismatchKernel, self).__init__(strings, **kwargs)
        self._m = kwargs.get("m", 1)
        if not self._m >= 0:
            raise ValueError("m must be >= 0")

    def _check_instance(self, instance, _, symbol, depth):
        i, offset, num_mismatches = instance
        string_symbol = self._entities[i][offset + depth]
        if string_symbol != symbol:
            num_mismatches += 1
        new_instance = (i, offset, num_mismatches)
        return new_instance, num_mismatches <= self._m

    def _to_instance(self, i, offset):
        # The third element is the number of mismatches survived so far
        return (i, offset, 0)

class ProfileKernel(_RecursivePrefixStringKernel):
    """Implementation of the PSSM-based string kernel [Kuang04]_.

    :param pssms: list of PSSMs of the form ``(residue, transp_score)''.
    :param k: kmer size, inclusive (default: 1).
    :param threshold: threshold mutation probability to count as a hit
        (default: 6.0)
    :param alphabet: list of valid symbols (default: AMINOACIDS).
    :param min_survivors_per_node: minimum number of instances to proceed lower
        within the prefix tree (default: 1)
    """
    def __init__(self, pssms, **kwargs):
        super(ProfileKernel, self).__init__(pssms, **kwargs)
        self._threshold = kwargs.get("threshold", 6.0)

    def _check_instance(self, instance, j, _, depth):
        i, offset, score = instance
        # Recall that more likely means (i) probability closer to 1, but (ii)
        # score (= neg log of probability) closer to 0. So to capture the most
        # likely mutant k-mers we want the total score not to increase too
        # much, and more specifically to remain below a threshold.
        score += self._entities[i][offset + depth][1][j]
        return (i, offset, score), score <= self._threshold

    def _to_instance(self, i, offset):
        # The third element is the total score of the mutations so far
        return (i, offset, 0.0)
