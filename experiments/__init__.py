# -*- coding: utf-8 -*-

import os
import numpy as np
from SPARQLWrapper import SPARQLWrapper, JSON
from sklearn.utils import check_random_state
from sklearn.svm import SVC
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc

import ocelot.ontology as O
from ocelot.scheduler import Scheduler
from ocelot import split_tr_ts

class _Experiment(object):
    """Base class for all experiments.

    Parameters
    ----------
    src : str
        Directory where the raw database data is held.
    dst : str
        Directory where the results will be held.
    endpoint : str
        URI of the SPARQL endpoint.
    default_graph : str
        URI of the default graph.
    force_update: bool, optional. (defaults to False)
        Whether to ignore cached results altogether.
    seed : int or np.random.RandomStream or None, optional. (defaults to None)
        RNG.
    """
    def __init__(self, src, dst, endpoint, default_graph, force_update = False,
                 seed = None, *args, **kwargs):
        self._rng = check_random_state(seed)
        try:
            os.mkdir(dst)
        except:
            # something will fail later on if dst does not exist
            pass
        if not (src and os.path.isdir(src)):
            raise ValueError("'{}' is not a valid directory".format(src))
        if not (dst and os.path.isdir(dst)):
            raise ValueError("'{}' is not a valid directory".format(dst))
        self.src, self.dst = src, dst
        if not endpoint:
            raise ValueError("no endpoint given.")
        if not default_graph:
            raise ValueError("no default graph given.")
        self.ep = SPARQLWrapper(endpoint)
        self.default_graph = default_graph
        if not self._check_graph(default_graph):
            raise ValueError("no graph '{}' at endpoint '{}'".format(default_graph, endpoint))
        self.force_update = force_update

        self._scheduler = Scheduler(self.dst)

    def _check_graph(self, graph):
        """Checks whether a graph exists."""
        ans = self.query(u"ASK WHERE {{ GRAPH <{default_graph}> {{ ?s ?p ?o }} }}")
        return ans[u"boolean"] == True

    def query(self, query):
        """Performs a query at the given endpoint."""
        prefixes = ""
        for shortcut, namespace in O.BINDINGS:
            prefixes += "PREFIX {}: <{}>\n".format(shortcut, unicode(namespace))
        query = prefixes + query.format(default_graph=self.default_graph)
        self.ep.setQuery(query)
        self.ep.setReturnFormat(JSON)
        return self.ep.query().convert()

    @staticmethod
    def _cast(d):
        if d[u"type"] in (u"uri", u"bnode", u"literal"):
            return d[u"value"]
        elif d[u"type"] == u"typed-literal":
            if d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#integer":
                return int(d[u"value"])
            elif d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#float":
                return float(d[u"value"])
            elif d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#double":
                return float(d[u"value"])
            elif d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#integer":
                return d[u"value"]
        raise NotImplementedError("can not cast '{}'".format(d.items()))

    def iterquery(self, query, n = None):
        ans = self.query(query)
        assert ans and len(ans) and "results" in ans
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self._cast(v) for k, v in bindings.iteritems() }
            if not n is None:
                assert len(bindings) == n, bindings
            yield bindings

    def evaluate_svm(self, folds, ys, kernel, C=1.0):
        results = []
        for k, (ts_indices, tr_indices) in enumerate(folds):
            model = SVC(C=C, kernel="precomputed", class_weight="auto")

            tr_ys, ts_ys = split_tr_ts(ys, tr_indices, ts_indices)
            tr_matrix, ts_matrix = split_tr_ts(kernel.compute(),
                                               tr_indices, ts_indices)
            model.fit(train_gram, train_ys)
            pr_ys = model.predict(test_gram)

            fpr, tpr, _ = roc_curve(ts_ys, pr_ys)
            pr, rc, f1, sup = precision_recall_fscore_support(ts_ys, pr_ys)

            results.append((sup, f1, pr, rc, auc(fpr, tpr)))
        return results

    def _compute_kernel(self, Kernel, *args, **kwargs):
        kernel = Kernel(*args, **kwargs)
        kernel.check_and_fixup(kwargs.get("tol", 1e-10))
        return kernel.compute(),

    def run(self):
        raise NotImplementedError()

from yip09 import *
from sgd import *
from cafa13 import *
from full import *
