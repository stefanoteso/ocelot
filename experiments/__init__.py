# -*- coding: utf-8 -*-

import os
import cPickle as pickle
import numpy as np
from SPARQLWrapper import SPARQLWrapper, JSON
from modshogun import CombinedKernel, CustomKernel
from modshogun import BinaryLabels, MKLClassification
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc

import ocelot.ontology as O
from ocelot.services import _cls

class _Experiment(object):
    """Base class for all experiments.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.

    .. todo::
        Add support for standard SVM.

    .. todo::
        Add support for SBR.
    """
    def __init__(self, src, dst, endpoint, default_graph, force_update = False,
                 *args, **kwargs):
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
            if n != None:
                assert len(bindings) == n, bindings
            yield bindings

    def _pickle(self, what, path):
        with open(os.path.join(self.dst, path), "wb") as fp:
            pickle.dump(what, fp)

    def _depickle(self, path):
        with open(os.path.join(self.dst, path), "rb") as fp:
            return pickle.load(fp)

    def _cached(self, f, relpath, *args, **kwargs):
        relpath += ".pickle"
        try:
            print _cls(self), ": loading '{}'".format(relpath)
            assert not self.force_update, "update forced by user"
            y = self._depickle(relpath)
        except Exception, e:
            print _cls(self), "|", e
            print _cls(self), ": computing '{}'".format(relpath)
            y = f(*args, **kwargs)
            self._pickle(y, relpath)
        return y

    def _cached_kernel(self, K, num, relpath, *args, **kwargs):
        path = os.path.join(self.dst, relpath)
        try:
            print _cls(self), ": loading '{}'".format(path)
            assert not self.force_update, "update forced by user"
            kernel = DummyKernel(path + ".txt", num = num, check_psd = True)
        except Exception, e:
            print _cls(self), "|", e
            print _cls(self), ": computing '{}'".format(path)
            kernel = K(*args, **kwargs)
            assert not kernel is None
            kernel.check_and_fixup(kwargs.get("tol", 1e-10))
            kernel.save(path + ".txt")
            kernel.draw(path + ".png")
        return kernel

    def _compute_kernels(self, infos, xs, xxs, tolerance = 1e-10):
        """Helper for computing kernels."""
        x_to_i = {x: i for i, x in enumerate(xs)}
        xx_indices = [(x_to_i[x1], x_to_i[x2]) for x1, x2 in xxs]

        for path, compute in infos:
            path = os.path.join(self.dst, path)
            try:
                print _cls(self), ": loading '{}'".format(path)
                kernel = DummyKernel(path + ".txt", num = len(xs),
                                     check_psd = True)
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                kernel = compute()
                assert not kernel is None
                kernel.check_and_fixup(tolerance)
                kernel.save(path + ".txt")
                kernel.draw(path + ".png")

            path += "-pairwise"
            try:
                print _cls(self), ": loading '{}".format(path)
                pairwise_kernel = DummyKernel(path + ".txt", num = len(xxs),
                                              check_psd = True)
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                pairwise_kernel = PairwiseKernel(xx_indices, kernel)
                pairwise_kernel.check_and_fixup(tolerance) 
                pairwise_kernel.save(path + ".txt")
                pairwise_kernel.draw(path + ".png")

            del kernel
            del pairwise_kernel

        kernels, pairwise_kernels = [], []
        for path, compute_kernel in infos:
            path = os.path.join(self.dst, path)
            print _cls(self), ": loading '{}".format(path)
            kernels.append(DummyKernel(path + ".txt", num = len(xs)))
            path += "-pairwise"
            print _cls(self), ": loading '{}".format(path)
            pairwise_kernels.append(DummyKernel(path + ".txt", num = len(xxs)))
        return kernels, pairwise_kernels

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
        combined_kernel = CombinedKernel()
        for matrix in matrices:
            combined_kernel.append_kernel(CustomKernel(matrix))
        return combined_kernel

    def _train_test_mkl(self, ys_tr, k_tr, ys_ts, k_ts, norm = 1.0, c = 1.0, class_costs = [0.5, 0.5]):
        """Run a single training/test round with MKL."""
        print _cls(self), ": running MKL (norm={} C={} costs={})" \
                            .format(norm, c, class_costs)

        # Create the model
        model = MKLClassification()
        model.set_C_mkl(c)
        model.set_mkl_norm(norm)
        model.set_C(*class_costs)

        # Learn the model
        model.set_kernel(k_tr)
        model.set_labels(BinaryLabels(ys_tr))
        model.train()

        # Infer
        k_ts.set_subkernel_weights(k_tr.get_subkernel_weights())
        model.set_kernel(k_ts)
        ys_pr = model.apply().get_labels()

        # Compute the results
        fpr, tpr, _ = roc_curve(ys_ts, ys_pr)
        pr_rc_f1_sup = precision_recall_fscore_support(ys_ts, ys_pr)
        return pr_rc_f1_sup + ( auc(fpr, tpr), )

    def _crossvalidate_mkl(self, folds, ys, kernels, **hyperparams):
        """Perform a k-fold coross-validation with MKL."""
        results = []
        for i, (ts_indices, tr_indices) in enumerate(folds):
            # Split the labels between training and test, and compute the
            # per-class costs on the training labels
            ys_tr, ys_ts = self._split_vector(ys, tr_indices, ts_indices)
            costs = self._get_class_costs(ys_tr)

            # Split the kernels between training and test
            matrices_tr, matrices_ts = [], []
            for kernel in kernels:
                matrix = kernel.compute()
                matrix_tr, matrix_ts = \
                    self._split_matrix(matrix, tr_indices, ts_indices)
                matrices_tr.append(matrix_tr)
                matrices_ts.append(matrix_ts)

            # Build the combined train/test kernels
            k_tr = self._combine_matrices(matrices_tr)
            k_ts = self._combine_matrices(matrices_ts)

            # Evaluate MKL
            result = self._train_test_mkl(ys_tr, k_tr, ys_ts, k_ts,
                                          class_costs = costs,
                                          **hyperparams)
            results.append(results)
        return results

    def _run_mkl(self, ys, kernels, folds):
        c_to_results = {}
        for c in (1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4):
            c_to_results[c] = self._crossvalidate_mkl(folds, ys, kernels,
                                                      c = c, norm = 1.0)
        from pprint import pprint
        pprint(c_to_results)

    def run(self):
        raise NotImplementedError

from yip09 import *
from sgd import *
from cafa13 import *
from full import *
