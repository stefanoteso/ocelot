# -*- coding: utf8 -*-

import numpy as np

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
    def is_psd(self, matrix):
        is_symmetric = (matrix.T == matrix).all()
        is_semi_psd = np.all(np.linalg.eigvalsh(matrix) >= 0)
        return is_symmetric and is_semi_psd
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

class DummyKernel(_Kernel):
    def __init__(self, arg, **kwargs):
        if isinstance(arg, str):
            matrix = np.loadtxt(arg)
        else:
            matrix = arg
        assert matrix.shape[0] == matrix.shape[1]
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

class SetKernel(_Kernel):
    """A generic set kernel.

    :param entities: list of dictionaries.
    """
    def _compute_all(self):
        num = len(self)
        matrix = np.zeros((num, num), dtype=np.float64)
        for i in xrange(num):
            set_i = self._entities[i]
            for j in xrange(i + 1):
                set_j = self._entities[j]
                if i == j and len(set_i) == 0 and len(set_j) == 0:
                    dp = 1.0
                else:
                    dp = float(len(set_i & set_j))
                matrix[i,j] = matrix[j,i] = dp
        return matrix

class ProfileKernel(_Kernel):
    """Profile-based string kernel.

    It acts as a wrapper around ``fastprofkernel``.

    *References*

    .. [Kuang04] Kuang et al., "Profile-based string kernels for remote
        homology detection and motif extraction", J. Bioinform. Comput. Biol.,
        2004.
    .. [Hamp13] Hamp et al., "Accelerating the Original Profile Kernel", PLoS
        One, 2013.
    .. [FPK] https://rostlab.org/owiki/index.php/Fastprofkernel 
    """
    def _compute_all(self, sequences):
        with open("sequences.fasta", "wt") as fp:
            for i, sequence in enumerate(sequences):
                fp.write(">{}\n{}\n".format(i, sequence))
        with open("identifiers.txt", "wt") as fp:
            for i in xrange(len(sequences)):
                fp.write(">{}\n".format(i))
        # XXX run fastprofkernel here
        profkernel = Binary("profkernel-core")
        args = [ "-i /usr/share/fastprofkernel/data/Amino.txt",
                 "-g /usr/share/fastprofkernel/data/Globals.txt",
                 "-o identifiers.txt",
                 "-L 4", "-Y 6.0", "-K",
                 "sequences.ascii-pssm" ]
        ret, out, err = profkernel.run(args)
        assert ret == 0

        return None

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
        out_i = 0
        for i, j in self._entities:
            out_j = 0
            for n, m in self._entities:
                kin, kim = submatrix[i,n], submatrix[i,m]
                kjn, kjm = submatrix[j,n], submatrix[j,m]
                matrix[out_i, out_j] = self.op(kin, kjm, kim, kjn)
                out_j += 1
            out_i += 1
        return matrix
