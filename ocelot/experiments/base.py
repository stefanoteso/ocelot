# -*- coding: utf-8 -*-

import numpy as np

import ocelot.ontology as O
from ocelot.services import _cls

class _Experiment(object):
    """Base class for all experiments.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    def __init__(self, src, dst, subtargets,
                 endpoint = u"http://127.0.0.1:8890/sparql",
                 default_graph = u"http://ocelot.disi.unitn.it/graph"):
        from SPARQLWrapper import SPARQLWrapper
        import os
        try:
            os.mkdir(dst)
        except:
            # something will fail later on if dst does not exist
            pass
        if not (src and os.path.isdir(src)):
            raise ValueError("'{}' is not a valid directory".format(src))
        if not (dst and os.path.isdir(dst)):
            raise ValueError("'{}' is not a valid directory".format(dst))
        self.src = src
        self.dst = dst
        self.subtargets = subtargets
        if not endpoint:
            raise ValueError("no endpoint given.")
        if not default_graph:
            raise ValueError("no default graph given.")
        self.ep = SPARQLWrapper(endpoint)
        self.default_graph = default_graph
        if not self._check_graph(default_graph):
            raise ValueError("no graph '{}' at endpoint '{}'".format(default_graph, endpoint))

    def _check_graph(self, graph):
        """Checks whether a graph exists."""
        ans = self.query(u"ASK WHERE {{ GRAPH <{default_graph}> {{ ?s ?p ?o }} }}")
        return ans[u"boolean"] == True

    def query(self, query):
        """Performs a query at the given endpoint."""
        from SPARQLWrapper import JSON
        prefixes = ""
        for shortcut, namespace in O.BINDINGS:
            prefixes += "PREFIX {}: <{}>\n".format(shortcut, unicode(namespace))
        query = prefixes + query.format(default_graph=self.default_graph)
        self.ep.setQuery(query)
        self.ep.setReturnFormat(JSON)
        return self.ep.query().convert()

    @staticmethod
    def cast(d):
        if d[u"type"] == u"uri":
            return d[u"value"]
        elif d[u"type"] == u"literal":
            return d[u"value"]
        elif d[u"type"] == u"typed-literal":
            if d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#integer":
                return int(d[u"value"])
            elif d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#integer":
                return float(d[u"value"])
            elif d[u"datatype"] == u"http://www.w3.org/2001/XMLSchema#integer":
                return d[u"value"]
        raise NotImplementedError("can not cast '{}'".format(d.items()))

    @staticmethod
    def _split_vector(ys, tr_indices, ts_indices):
        return ys[np.ix_(tr_indices)], ys[np.ix_(ts_indices)]

    @staticmethod
    def _split_matrix(k, tr_indices, ts_indices, mode = "shogun"):
        k_tr = k[np.ix_(tr_indices, tr_indices)]
        if mode == "shogun":
            k_ts = k[np.ix_(tr_indices, ts_indices)]
        elif mode == "sklearn":
            k_ts = k[np.ix_(ts_indices, tr_indices)]
        else:
            raise ValueError("invalid mode '{}'".format(mode))
        return k_tr, k_ts

    @staticmethod
    def _get_class_costs(ys):
        num_pos = float(len(filter(lambda y: y >= 0, ys)))
        return (len(ys) - num_pos) / len(ys), num_pos / len(ys)

    @staticmethod
    def _combine_matrices(matrices):
        from modshogun import CombinedKernel, CustomKernel
        combined_kernel = CombinedKernel()
        for matrix in matrices:
            combined_kernel.append_kernel(CustomKernel(matrix))
        return combined_kernel

    def _crossvalidate_mkl(self, ys, kernels, folds, mkl_c = 1.0, mkl_norm = 1):
        from modshogun import BinaryLabels, MKLClassification
        from sklearn.metrics import precision_recall_fscore_support

        results = []
        for i, (ts_indices, tr_indices) in enumerate(folds):
            print _cls(self), ": fold {}/{}, preparing; C={} norm={}".format(i, len(folds), mkl_c, mkl_norm)

            # Split the labels between training and test, and compute the
            # per-class costs on the training labels
            ys_tr, ys_ts = self._split_vector(ys, tr_indices, ts_indices)
            cost_pos, cost_neg = self._get_class_costs(ys_tr)

            # Split the kernels between training and test
            matrices_tr, matrices_ts = [], []
            for kernel in kernels:
                matrix = kernel.compute()
                matrix_tr, matrix_ts = self._split_matrix(matrix, tr_indices, ts_indices)
                matrices_tr.append(matrix_tr)
                matrices_ts.append(matrix_ts)
            combined_kernel_tr = self._combine_matrices(matrices_tr)
            combined_kernel_ts = self._combine_matrices(matrices_ts)

            # Create the model
            model = MKLClassification()
            model.set_C_mkl(mkl_c)
            model.set_mkl_norm(mkl_norm)
            model.set_C(cost_pos, cost_neg)

            # Learn
            print _cls(self), ": fold {}/{}, learning (class costs = +{} -{})".format(i, len(folds), cost_pos, cost_neg)
            model.set_kernel(combined_kernel_tr)
            model.set_labels(BinaryLabels(ys_tr))
            model.train()
            beta = combined_kernel_tr.get_subkernel_weights()

            # Infer
            print _cls(self), ": fold {}/{}, predicting".format(i, len(folds))
            model.set_kernel(combined_kernel_ts)
            combined_kernel_ts.set_subkernel_weights(beta)
            ys_pr = model.apply().get_labels()

            results.append(precision_recall_fscore_support(ys_ts, ys_pr))

        return results

    def _run_mkl(self, ys, kernels, folds):
        c_to_results = {}
        for c in (1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4):
            c_to_results[c] = self._crossvalidate_mkl(ys, kernels, folds, mkl_c = c)
        from pprint import pprint
        pprint(c_to_results)

    def _crossvalidate_sbr_with_mkl(self, ys, ks, folds):
        pass

    def run(self):
        raise NotImplementedError

