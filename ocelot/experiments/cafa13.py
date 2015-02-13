# -*- coding: utf-8 -*-

from .base import _Experiment

class CAFA13Experiment(_Experiment):
    """New experiment based on the CAFA13 dataset."""
    def __init__(self, *args, **kwargs):
        super(YeastExperiment, self).__init__(*args, **kwargs)
    def run(self):
        raise NotImplementedError()
