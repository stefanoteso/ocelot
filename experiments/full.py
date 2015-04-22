# -*- coding: utf-8 -*-

from . import _Experiment

class AllSpeciesExperiment(_Experiment):
    """New experiment based of species-agnostic BioGRID and iPfam."""
    def __init__(self, *args, **kwargs):
        super(AllSpeciesExperiment, self).__init__(*args, **kwargs)
    def run(self):
        raise NotImplementedError()
