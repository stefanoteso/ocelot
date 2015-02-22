# -*- coding: utf8 -*-

import numpy as np
from ocelot.services import _cls

class _Kernel(object):
    """Base kernel class.

    :param entities: list of arbitrary data entries.
    :param do_normalize: normalize the Gram matrix.
    """
    def __init__(self, entities, **kwargs):
        self._entities = entities
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
    def draw(self, path):
        """Draws the kernel matrix to a file."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm
        except:
            print "matplotlib is required; can not draw"
            return
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        matrix = self.compute()
        ax.matshow(matrix, interpolation = "nearest", cmap = cm.OrRd)
        fig.savefig(path)
    def save(self, path):
        np.savetxt(path, self.compute())
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

class DummyKernel(_Kernel):
    def __init__(self, arg, **kwargs):
        num = kwargs.get("num")
        if isinstance(arg, str):
            matrix = np.loadtxt(arg)
        else:
            matrix = arg
        assert matrix.shape[0] == matrix.shape[1]
        if num != None:
            assert matrix.shape == (num, num)
        super(DummyKernel, self).__init__(range(matrix.shape[0]),
                                          **kwargs)
        self._matrix = matrix
        if self._do_normalize:
            self._matrix = self._normalize(self._matrix)

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

class SpectrumKernel(_Kernel):
    """The spectrum kernel.

    NOTE: this is a very dumb, slow implementation. Pro implementation use
    tries, kids!

    :param strings: a list of strings.
    :param kmin: minimum kmer size (inclusive).
    :param kmax: maximum kmer size (exclusive).
    :param binary: check for matches rather than counting them.

    *References*

    [1] Leslie et al. "The spectrum kernel: A string kernel for SVM protein
        classification.", 2002.
    """
    def __init__(self, strings, **kwargs):
        super(SpectrumKernel, self).__init__(strings, **kwargs)
        self.kmin = kwargs.get("kmin", 1)
        self.kmax = kwargs.get("kmax", self.kmin + 1)
        if self.kmin < 1 or self.kmin > self.kmax:
            raise ValueError, "invalid kmin or kmax, '{}' '{}'".format(self.kmin, self.kmax)
        self.binary = kwargs.get("binary", False)
    def _count_substrings(self, string, k):
        n = len(string)
        if k > n:
            return {}
        counts = {}
        for i in xrange(0, n - k + 1):
            substring = string[i:i+k]
            if not substring in counts:
                counts[substring] = 1
            else:
                counts[substring] += 1
        return counts
    def _compute_all(self):
        n = len(self)
        matrix = np.zeros((n, n))
        for k in xrange(self.kmin, self.kmax):
            # First we collect all substrings of all strings
            print "counting k=", k
            counts = {}
            for string in self._entities:
                counts[string] = self._count_substrings(string, k)
            print "computing k=", k
            # Now we compare them
            for i in xrange(n):
                counts_i = counts[self._entities[i]]
                for j in xrange(i, n):
                    counts_j = counts[self._entities[j]]
                    common_substrings = set(counts_i.keys()) & set(counts_j.keys())
                    dp = 0.0
                    if self.binary:
                        dp += len(common_substrings)
                    else:
                        for substring in common_substrings:
                            dp += counts_i[substring] * counts_j[substring]
                    matrix[i,j] = matrix[j,i] = dp
        return matrix

class MismatchKernel(_Kernel):
    """Mismatch string kernel."""
    def __init__(self, strings, k = 3, m = 2, **kwargs):
        super(MismatchKernel, self).__init__(strings, **kwargs)
    def _get_neighborhood(self, alpha):
        raise NotImplementedError 
    def _compute_all(self):
        all_counts = []
        for string in self._entities:
            counts = {}
            for i in xrange(len(string) - k + 1):
                alpha = string[i:i+k]
                for beta in self._get_neighborhood(alpha):
                    counts[beta] += 1
            all_counts.append(counts)
        return SparseLinearKernel(all_counts).compute()

# Taken from fastprofkernel
_BACKGROUND_FREQ = {
    "A": 0.0799912015849807,
    "C": 0.0171846021407367,
    "D": 0.0578891399707563,
    "E": 0.0638169929675978,
    "F": 0.0396348024787477,
    "G": 0.0760659374742852,
    "H": 0.0223465499452473,
    "I": 0.0550905793661343,
    "K": 0.060458245507428,
    "L": 0.0866897071203864,
    "M": 0.0215379186368154,
    "N": 0.044293531582512,
    "P": 0.0465746314476874,
    "Q": 0.0380578923048682,
    "R": 0.0484482507611578,
    "S": 0.0630028230885602,
    "T": 0.0580394726014824,
    "V": 0.0700241481678408,
    "W": 0.0144991866213453,
    "Y": 0.03635438623143,
}

class ProfileKernel(_Kernel):
    """Profile-based string kernel [Kuang04]_.

    Please note that there are faster implementations around, namely the
    ``fastprofkernel`` package [fastprofkernel]_ by RostLab.

    :param pssms: list of PSSM matrices.
    :param k: length of the k-mers.
    :param threshold: threshold mutation probability to count as a hit.

    *References*

    .. [Kuang04] Kuang et al., "Profile-based string kernels for remote
        homology detection and motif extraction", J. Bioinform. Comput. Biol.,
        2004.

    .. [fastprofkernel] `<https://rostlab.org/owiki/index.php/Fastprofkernel>`
    """ 
    def __init__(self, pssms, k = 4, threshold = 6.0, **kwargs):
        super(ProfileKernel, self).__init__(pssms, **kwargs)
        self.k, self.m = k, m
        self.threshold = threshold
    def _compute_all(self):
        all_counts = []
        for string in self._entities:
            counts = {}
            for i in xrange(len(string) - k + 1):
                kmer = string[i:i+k]
                for mmer in self._get_neighborhood(kmer):
                    counts[mmer] += 1
            all_counts.append(counts)
        return SparseLinearKernel(all_counts).compute()

class DiffusionKernel(_Kernel):
    """The diffusion kernel between graph nodes.

    *References*

    .. [DK] Kondor and Lafferty, "Diffusion Kernels on Graphs and Other
        Discrete Input Spaces", 2002.
    """
    def __init__(self, adjmatrix, beta = 1.0, **kwargs):
        n, m = adjmatrix.shape
        assert n == m
        self._adjmatrix = adjmatrix
        self._beta = beta
        super(DiffusionKernel, self).__init__(range(n), **kwargs)
    def _compute_all(self):
        import scipy.linalg as la
        # XXX this does not work for weighted adjacency matrices
        a = self._adjmatrix
        d = np.identity(a.shape[0]) * np.sum(a, axis = 1)
        e = la.expm(self._beta * (a - d))

        # being approximate, matrix exponentiation can give slightly asymmetric
        # results; here we fix them if they are small enough.
        asymmetry = np.linalg.norm(e.T - e, ord = "fro")
        if asymmetry > 0.0:
            assert asymmetry < 1e-10, "the gods are playful today"
            e = 0.5 * (e.T + e)
        return e

class EmpiricalKernelMap(_Kernel):
    """Empirical kernel map over an underlying kernel[1].

    *References*

    [1] Scholkopf et al., "Input Space Versus Feature Space in Kernel-Based
        Methods", 1999.
    """
    def __init__(self, *args, **kwargs):
        import ocelot.features
        super(EmpiricalKernelMap, self).__init__(*args, **kwargs)
        self.features = features.EmpiricalFeatures(args[0], range(len(self)))
    def _compute_all(self):
        LinearKernel(self.features, range(len(self)))

class PairwiseKernel(_Kernel):
    """The pairwise kernel.

    It computes a kernel between pairs of elements as follows:

    .. math::
        k[(i,j),(n,m)] := k[i,n]k[j,m] + k[i,m]k[j,n]

    or:

    .. math::
        k[(i,j),(n,m)] := k[i,n] + k[j,m] + k[i,m] + k[j,n]

    :param indices: indices of the elements to combine.
    :param subkernel: a `Kernel` instance.
    :param op: either ``"tensor"`` or ``"sum"``.
    """
    def __init__(self, pairs, subkernel, op = "tensor", **kwargs):
        super(PairwiseKernel, self).__init__(pairs, **kwargs)
        self.subkernel = subkernel
        if op == "tensor":
            self.op = lambda kin, kjm, kim, kjn: kin * kjm + kim * kjn
        elif op == "sum":
            self.op = lambda kin, kjm, kim, kjn: kin + kjm + kim + kjn
        else:
            raise ValueError, "invalid op '{}'".format(op)
    def _compute_all(self):
        submatrix = self.subkernel.compute()
        matrix = np.zeros((len(self), len(self)))
        for out_i, (i, j) in enumerate(self._entities):
            for out_j, (n, m) in enumerate(self._entities):
                kin, kim = submatrix[i,n], submatrix[i,m]
                kjn, kjm = submatrix[j,n], submatrix[j,m]
                matrix[out_i, out_j] = self.op(kin, kjm, kim, kjn)
        return matrix
