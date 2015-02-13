# -*- coding: utf-8 -*-

from .base import _Experiment

class YeastExperiment(_Experiment):
    """New experiment based on SGD and iPfam."""
    def __init__(self, *args, **kwargs):
        super(YeastExperiment, self).__init__(*args, **kwargs)
    def run(self):
        raise NotImplementedError()
    # figure out the PPI, GO, SEQ
    # figure out the DDI-RRI (only use instance information)
    # match the PPI and DDI-RRI
    # compute composition features
    # **compute complexity features
    # **compute conservation (profile) features
    # **compute secondary features
    # compute subloc features (binary vec)
    # compute COG features (binary vec)
    # compute cell-cycle gene expression features (correlations)
    # compute environment-response gene expression features (correlations)
    # read in the Y2H raw data
    # read in the TAP-MS raw data
    # read in the Yip folds
    # write the SBR output (individuals, kernels, predicates, rules)
    # TODO: filter away short proteins
    # TODO: redundancy-reduce the protein dataset, map proteins to their
    #       representatives
    # TODO: form the CV folds
    # TODO: put all duplicates back into the training set
