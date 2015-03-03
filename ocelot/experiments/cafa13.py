# -*- coding: utf-8 -*-

from .base import _Experiment

class CAFA13Experiment(_Experiment):
    """New experiment based on the CAFA13 dataset.

    This experiment only targets protein function prediction.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    def __init__(self, *args, **kwargs):
        super(CAFA13Experiment, self).__init__(*args, **kwargs)

    def run(self):
        raise NotImplementedError()
