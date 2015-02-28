# -*- coding: utf8 -*-

import multiprocessing as mp
import numpy as np
from scipy import sparse
from ocelot.services import AMINOACIDS

from .base import _Kernel

MAX_SURVIVORS_PER_FLUSH = 65536

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
    :param num_processes: number of processes to distribute the computation on,
        0 means *use all CPUs* (default: 1)

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
        if not len(self._alphabet) >= 2:
            raise ValueError("alphabet must have at least two symbols")
        self._min_survivors_per_node = kwargs.get("min_survivors_per_node", 1)
        if not self._min_survivors_per_node >= 1:
            raise ValueError("min_survivors_per_node must be >= 1")
        self._num_processes = kwargs.get("num_processes", 1)
        if self._num_processes == 0:
            self._num_processes = mp.cpu_count()
        if not self._num_processes >= 1:
            raise ValueError("num_processes must be >= 1")
        # Used for the delayed-update optimization
        self._survivors_rows, self._survivors_cols, self._survivors_counts = \
            [], [], []
        # Precompute a couple things
        self._num_nodes = len(self._alphabet)**(self._k + 1) - 1
        self._index_symbol_pairs = list(enumerate(self._alphabet))

    def _flush_survivors(self):
        assert len(self._survivors_rows) == len(self._survivors_counts) and \
               len(self._survivors_cols) == len(self._survivors_counts)
        v = sparse.coo_matrix((self._survivors_counts,
                               (self._survivors_rows, self._survivors_cols)),
                              shape = (len(self._entities), self._num_nodes)).tocsr()
        self._matrix += v * v.T
        self._survivors_rows, self._survivors_cols, self._survivors_counts = \
            [], [], []

    def _add_survivors(self, instances):
        """Called upon visiting a leaf of the prefix tree. Appends the list of
        survivors with the surviving instances."""
        for i in xrange(len(self._entities)):
            count = sum(1 for instance in instances if instance[0] == i)
            if count > 0:
                self._survivors_rows.append(i)
                self._survivors_cols.append(self._node_counter)
                self._survivors_counts.append(count)
        if len(self._survivors_counts) > MAX_SURVIVORS_PER_FLUSH:
            self._flush_survivors()

    def _recurse(self, instances, depth):
        """Depth-first traversal of the prefix-tree."""
        assert 0 <= depth <= self._k
        self._node_counter += 1
        if depth == self._k:
            self._add_survivors(instances)
        else:
            for s, symbol in self._index_symbol_pairs:
                # Filter out all instances that do not pass the check at the
                # current node 
                survivors = []
                for instance in instances:
                    new_instance, check = \
                        self._check_instance(instance, s, symbol, depth)
                    if check:
                        survivors.append(new_instance)
                # If there are no instances left, there is nothing else to do
                if len(survivors) < self._min_survivors_per_node:
                    continue
                # Process the child nodes
                self._recurse(survivors, depth + 1)
                # Get rid of the list of survivors
                del survivors

    def _to_instance(self, i, offset):
        """Given an element index and an offest, returns an instance (i.e.
        the k-mer indices optionally associated to additional information).

        **Note**: the first element of an instance **must** be the index of
        the pssm it refers to! The ``_add_survivors`` method relies on it.
        """
        return (i, offset)

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

        if self._num_processes > 1:
            pool = mp.Pool(self._num_processes)
            pool.map(self._recurse, self._instances)
        else:
            self._node_counter = 0
            self._recurse(self._instances, 0)
            self._flush_survivors()
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
    """The profile kernel for strings [Kuang04]_.

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

class _TestRecursivePrefixStringKernel(object):
    def test_k_too_small(self):
        import pytest
        strings = ("TEST", "TEST")
        with pytest.raises(ValueError):
            kernel = _RecursivePrefixStringKernel(strings, k = 0)
    def test_k_too_large(self):
        import pytest
        for k in xrange(10):
            string = "A" * k
            strings = (string, string)
            with pytest.raises(ValueError):
                kernel = _RecursivePrefixStringKernel(strings, k = k + 1)
                instances = kernel._get_instances()
    def test_alphabet(self):
        import pytest
        strings = ("TEST", "TEST")
        with pytest.raises(ValueError):
            kernel = _RecursivePrefixStringKernel(strings, k = 0, alphabet = "")
        with pytest.raises(ValueError):
            kernel = _RecursivePrefixStringKernel(strings, k = 0, alphabet = "?")
    def test_num_instances(self):
        from ocelot.services import AMINOACIDS
        STRING = "".join(AMINOACIDS)
        for num_replicas in xrange(2, 5):
            strings = [STRING] * num_replicas
            for k in xrange(1, len(STRING)):
                kernel = _RecursivePrefixStringKernel(strings, k = k)
                instances = kernel._get_instances()
                assert len(instances) == num_replicas * (len(STRING) - k + 1)

class _TestSpectrumKernel(object):
    def test_result(self):
        STRINGSETS = (
            (("A", "A"), 1, np.array([[1, 1], [1, 1]])),
            (("A", "Y"), 1, np.array([[1, 0], [0, 1]])),
            (("AA", "AA"), 1, np.array([[4, 4], [4, 4]])),
            (("AA", "YY"), 1, np.array([[4, 0], [0, 4]])),
            (("AA", "AA"), 2, np.array([[1, 1], [1, 1]])),
            (("AA", "YY"), 2, np.array([[1, 0], [0, 1]])),
            (("AY", "YA"), 1, np.array([[2, 2], [2, 2]])),
            (("AY", "AA"), 1, np.array([[2, 2], [2, 4]])),
        )
        for strings, k, expected in STRINGSETS:
            kernel = SpectrumKernel(strings, k = k, do_normalize = False)
            output = kernel.compute()
            assert (output == expected).all()

class _TestMismatchKernel(object):
    def test_foo(self):
        # WRITEME
        pass

class _TestProfileKernel(object):
    def test_foo(self):
        # WRITEME
        pass

