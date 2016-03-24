# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.stats import entropy
from services import *

# TODO add B (Asparagine or Aspartic Acid), Z (Glutamine or Glutamic Acid), J
# (Leucine or Isoleucine), U (Selenocystein), O (Pyrrolisine)
AA_INFO = pd.read_csv("data/aainfo.tsv", sep="	").set_index(u"aa")
"""Amino acid information, taken from [WikipediaAAa]_, [WikipediaAAb]_,
[Swart04]_."""

AA1_TO_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}
AA3_TO_AA1 = {aa3: aa1 for aa1, aa3 in AA1_TO_AA3.iteritems()}

AMINOACIDS = list(AA_INFO.index)
AMINOACIDS3 = [AA1_TO_AA3[aa] for aa in AMINOACIDS]

AA_INFO.loc["?"] = AA_INFO.mean()

def get_composition_phi(sequence, what=None):
    """Composition features based on ``AA_INFO``.

    Parameters
    ----------
    sequence : str
        The amino acid sequences.
    what : list, optional. (defaults to ``None``)
        Target composition information. ``None`` means all of them.

    Returns
    -------
    phi : np.ndarray of shape (nres, nwhat)
        Array with one residue per row, one target feature per column.
    """
    if what is None:
        what = list(AA_INFO.columns)
    aas = [aa if aa in AA_INFO.index else "?" for aa in sequence.upper()]
    print aas
    return AA_INFO.loc[aas][what].values

def get_low_complexity_phi(sequence, windows=None, epsilon=1e-10):
    """Low-complexity region features.

    Low-complexity regions are correlated with packing density and
    hydrophobicity and with the interface propensity of residues.
    Some useful references are [Wootton94]_, [Liao05]_, [Coletta10]_.

    Measured by computing the Shannon entropy over one or more windows of
    residues.

    Parameters
    ----------
    sequence :
        The amino acid sequence.
    windows : list, optional. (defaults to [8, ..., 128])
        WRITEME
    epsilon : float, optional. (defaults to 1e-10)
        WRITEME
    """
    if windows is None:
        windows = [8, 16, 32, 64, 128]
    phi = []
    for window in windows:
        window_phi = []

        padding = "?" * window
        s = padding + sequence + padding
        for i in range(len(sequence)):
            i_min = len(padding) + i - window
            i_max = len(padding) + i + window + 1
            w = s[i_min:i_max]

            counts = [len(filter(lambda res: res == aa, w)) + epsilon
                      for aa in AMINOACIDS]

            h = entropy(np.array(counts, dtype=np.float32) / sum(counts))
            if np.abs(h) < epsilon:
                h = 0.0
            window_phi.append(h)

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

class EmpiricalFeatures(object):
    """The Empirical Kernel Map [Scholkopf99]_.

    Given a fixed set of patterns :math:`(z_1, \\ldots, z_m)` and an input
    pattern :math:`x`, the empirical kernel map is defined as:

    .. math::

        \\varphi(x) := (k(x,z_1), ..., k(x,z_m))

    for any given sub-kernel :math:`k`.

    NOTE: patterns and indices refer to elements in the sub-kernel.
    """
    def __init__(self, subkernel, patterns):
        self.subkernel = subkernel
        self.patterns = patterns
    def compute(self, indices):
        return self.subkernel.compute()[indices, self.patterns]
