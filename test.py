#!/usr/bin/env python2
# -*- coding: utf8 -*-

import unittest
import numpy as np

from ocelot import SpectrumKernel, MismatchKernel, AMINOACIDS

_AAS = "".join(AMINOACIDS)

_STRINGSETS = (
    ("AA", "YY"),
    ("AA", "AA", "YY", "YY"),
    ("AAA", "AAA", "YYY", "YYY"),
    ("AA", "AY", "YA", "YY"),
    (_AAS[:5], _AAS[:5]),
)

class SpectrumKernelTest(unittest.TestCase):
    def test_kmers(self):
        kernel = SpectrumKernel(_STRINGSETS[0], k = 1)
        kernel.compute()
    def test_results(self):
        for stringset in _STRINGSETS:
            max_len = max(map(len, stringset))
            for k in xrange(1, max_len + 1):
                kernel = SpectrumKernel(stringset, k = k, do_normalize = False)
                matrix = kernel.compute()
                self.assertTrue(matrix.shape == (len(stringset), len(stringset)))

class MismatchKernelTest(unittest.TestCase):
    def test_kmers(self):
        kernel = MismatchKernel(_STRINGSETS[0], k = 1, m = 1)
        kernel.compute()
    def test_results(self):
        for stringset in _STRINGSETS:
            max_len = max(map(len, stringset))
            for k in xrange(1, max_len + 1):
                kernel = MismatchKernel(stringset, k = k, m = 1, do_normalize = False)
                matrix = kernel.compute()
                self.assertTrue(matrix.shape == (len(stringset), len(stringset)))

def _pssm1(i):
    assert 0 <= i <= 19
    pssm = [ 0.0 for _ in xrange(20) ]
    pssm[i] = 1.0
    return pssm

_PROFILESETS = (
    (("A", _pssm1(0)),
     ("A", _pssm1(0))),
    (("Y", _pssm1(5)),
     ("Y", _pssm1(5))),
)

class ProfileKernelTest(unittest.TestCase):
    def test_kmers(self):
        kernel = MismatchKernel(_PROFILESETS[0], k = 1, )
        kernel.compute()
    def test_results(self):
        for stringset in _PROFILESETS:
            max_len = max(map(len, stringset))
            for k in xrange(1, max_len + 1):
                pass

if __name__ == "__main__":
    unittest.main() 
