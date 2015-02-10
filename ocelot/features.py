# -*- coding: utf-8 -*-

import numpy as np
from collections import namedtuple
from services import *

# Amino acid one- and three-letter codes
_OLC_TO_TLC = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}
_TLC_TO_OLC = { v: k for k, v in _OLC_TO_TLC.items() }
_OLCS = _OLC_TO_TLC.keys()
_TLCS = _OLC_TO_TLC.values()

# Amino acid composition
# XXX add B (Asparagine or Aspartic Acid), Z (Glutamine or Glutamic Acid), J
#     (Leucine or Isoleucine), U (Selenocystein), O (Pyrrolisine)
_Info = namedtuple("_Info", ("polarity", "charge", "hydropathy",
                             "hydrophobicity", "avgmass1", "avgmass2",
                             "volume", "isoelectric_point", "pk1", "pk2",
                             "aromaticity", "polarizability"))
_COMPOSITION = {
    "A": _Info(-1.0,  0.0,  1.8,  1.0,   89.09,  71.08,  67.0,  6.01, 2.35,  9.87,  0.0,  50.16),
    "C": _Info(-1.0,  0.0,  2.5,  1.0,  121.15, 103.14,  86.0,  5.05, 1.92, 10.70,  0.0,  74.99),
    "D": _Info( 1.0, -1.0, -3.5, -1.0,  133.10, 115.09,  91.0,  2.85, 1.99,  9.90,  0.0,  69.09),
    "E": _Info( 1.0, -1.0, -3.5, -1.0,  147.13, 129.12, 109.0,  3.15, 2.10,  9.47,  0.0,  84.67),
    "F": _Info(-1.0,  0.0,  2.8,  1.0,  165.19, 147.18, 135.0,  5.49, 2.20,  9.31, -1.0, 121.43),
    "G": _Info(-1.0,  0.0, -0.4,  1.0,   75.07,  57.05,  48.0,  6.06, 2.35,  9.78,  0.0,  36.66),
    "H": _Info( 1.0,  0.1, -3.2, -1.0,  155.16, 137.14, 118.0,  7.60, 1.80,  9.33, -1.0,  99.48),
    "I": _Info(-1.0,  0.0,  4.5,  1.0,  131.17, 113.16, 124.0,  6.05, 2.32,  9.76,  1.0,  91.21),
    "K": _Info( 1.0,  1.0, -3.9, -1.0,  146.19, 128.17, 135.0,  9.60, 2.16,  9.06,  0.0, 101.73),
    "L": _Info(-1.0,  0.0,  3.8,  1.0,  131.17, 113.16, 124.0,  6.01, 2.33,  9.74,  1.0,  91.60),
    "M": _Info(-1.0,  0.0,  1.9,  1.0,  149.21, 131.20, 124.0,  5.74, 2.13,  9.28,  0.0, 102.31),
    "N": _Info( 1.0,  0.0, -3.5, -1.0,  132.12, 114.10,  96.0,  5.41, 2.14,  8.72,  0.0,  73.15),
    "P": _Info(-1.0,  0.0, -1.6,  1.0,  115.13,  97.12,  90.0,  6.40, 1.95, 10.64,  0.0,  73.47),
    "Q": _Info( 1.0,  0.0, -3.5, -1.0,  146.15, 128.13, 114.0,  5.65, 2.17,  9.13,  0.0,  88.79),
    "R": _Info( 1.0,  1.0, -4.5, -1.0,  174.20, 156.19, 148.0, 10.76, 1.82,  8.99,  0.0, 114.81),
    "S": _Info( 1.0,  0.0, -0.8, -1.0,  105.09,  87.08,  73.0,  5.68, 2.19,  9.21,  0.0,  53.82),
    "T": _Info(-1.0,  0.0, -0.7, -1.0,  119.12, 101.11,  93.0,  5.60, 2.09,  9.10,  0.0, 153.06),
    "V": _Info(-1.0,  0.0,  4.2,  1.0,  117.15,  99.13, 105.0,  6.00, 2.39,  9.74,  1.0,  76.09),
    "W": _Info( 1.0,  0.0, -0.9, -1.0,  204.23, 186.21, 163.0,  5.89, 2.46,  9.41, -1.0, 153.06),
    "Y": _Info( 1.0,  0.0, -1.3, -1.0,  181.19, 163.18, 141.0,  2.46, 2.20,  9.21, -1.0, 126.19),
}
_COMPOSITION_TARGETS = _COMPOSITION.values()[0]._fields

averages = []
for target in _COMPOSITION_TARGETS:
    averages.append(sum([ _COMPOSITION[aa]._asdict()[target]
                          for aa in _OLCS ]) / len(_COMPOSITION))
_COMPOSITION["?"] = _Info(*averages)

class CompositionFeatures(object):
    """Per-residue composition features for a single peptide sequence.

    :param targets: list of target composition features (strings).

    *References*

    .. [Compos1] https://en.wikipedia.org/wiki/Amino_acid
    .. [Compos2] https://en.wikipedia.org/wiki/Proteinogenic_amino_acid
    .. [Swart04] Swart, Snijders and van Duijnen, "Polarizabilities of amino
        acid residues", 2004.
    """
    def __init__(self, targets = None):
        if targets == None:
            targets = _COMPOSITION_TARGETS
        self._targets = targets
    def compute(self, sequence):
        """Computes the per-residue features.

        :param sequence: a polypeptide sequence of `N` residues.
        :returns: a `N*K` matrix, with `K` target features per row.
        """
        if not len(sequence):
            raise ValueError("empty sequence")
        phis = []
        for residue in sequence.upper():
            phi = [ _COMPOSITION[residue]._asdict()[target]
                    for target in self._targets]
            phis.append(phi)
        return np.array(phis)

class ComplexityFeatures(object):
    """Low-complexity region features.

    Low-complexity regions are correlated with packing density and
    hydrophobicity and with the interface propensity of residues.

    Measured by computing the Shannon entropy over one or more windows of
    residues.

    XXX implement priors (pseudocounts).

    *References*

    .. [Wootton94] Wootton J, "Sequences with unusual amino acid compositions",
        1994.
    .. [Liao05] Liao, Yeh, Chiang, Jernigan, Lustig, "Protein sequence entropy
        is closely related to packing density and hydrophobicity", 2005.
    .. [Coletta10] Coletta, Pinney, Weiss, Solís, Marsh, Pettifer, Attwood,
        "Low-complexity regions within protein sequences have
        position-dependent roles", 2010.
    """
    def compute(self, sequence, window_sizes = [8,16,32,64,128], epsilon = 1e-10):
        from scipy.stats import entropy
        import math
        phi = []
        for window_size in window_sizes:
            window_phi = []
            padding = "?" * window_size
            sequence_ = padding + sequence + padding
            for i in xrange(len(sequence)):
                i_min = len(padding) + i - window_size
                i_max = len(padding) + i + window_size + 1
                window = sequence_[i_min:i_max]
                counts = [ float(len(filter(lambda res: res == aa, window))) + epsilon
                           for aa in _OLCS ]
                window_h = entropy(np.array(counts) / np.sum(counts))
                if math.fabs(window_h) < epsilon:
                    window_h = 0.0
                window_phi.append(window_h)
            phi.append(window_phi)
        return np.array(phi).T

class SequenceFeatures(object):
    """Per-sequence features that wrap per-residue features.

    :param Features: wrapped per-residue feature class.

    Additional parameters are forwarded to the wrapped class.
    """
    def __init__(self, Features, *args, **kwargs):
        self.features = Features(*args, **kwargs)
    def compute(self, sequences):
        """Computes the per-sequence features.

        :param sequences: `M` polypeptide sequence.
        :returns: a `M*K` matrix, with `K` target features per row.
        """
        return np.array([ np.mean(self.features.compute(sequence), axis = 0)
                          for sequence in sequences ])

class SupFamFeatures(object):
    """Superfamily domain features.

    XXX requires UNIPROT IDs!

    Adapted from an implementation by Luca Masera.
    """
    def compute(self, sequences):
        pass

class EmpiricalFeatures(object):
    """The Empirical Kernel Map[1]

    Given a set of fixed patterns `(z1, ..., zm)` and an input pattern `x`, the
    empirical kernel map is defined as:

    .. math::

        phi(x) := (k(x,z1), ..., k(x,zm))

    for any given sub-kernel k.

    NOTE: patterns and indices refer to elements in the sub-kernel.

    *References*

    [1] Scholkopf et al., "Input Space Versus Feature Space in Kernel-Based
        Methods", 1999.
    """
    def __init__(self, subkernel, patterns):
        self.subkernel = subkernel
        self.patterns = patterns
    def compute(self, indices):
        return self.subkernel.compute()[indices, self.patterns]
