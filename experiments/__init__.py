# -*- coding: utf-8 -*-

import os
import cPickle as pickle
import numpy as np
from SPARQLWrapper import SPARQLWrapper, JSON
from modshogun import CombinedKernel, CustomKernel
from modshogun import BinaryLabels, MKLClassification
from sklearn.svm import SVC
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc
from collections import namedtuple

import ocelot.ontology as O
from ocelot.services import _cls

Stage = namedtuple("Stage", ["f", "inputs", "outputs"])

class _Experiment(object):
    """Base class for all experiments.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    :param force_update: whether to discard the cache contents (default: ``False``).

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
            if not n is None:
                assert len(bindings) == n, bindings
            yield bindings

    def _pickle(self, what, path):
        with open(os.path.join(self.dst, path), "wb") as fp:
            pickle.dump(what, fp)

    def _depickle(self, path):
        with open(os.path.join(self.dst, path), "rb") as fp:
            return pickle.load(fp)

#    def _cached(self, f, relpath, *args, **kwargs):
#        relpath += ".pickle"
#        try:
#            assert not self.force_update, "forced update"
#            print _cls(self), ": loading '{}'".format(relpath)
#            y = self._depickle(relpath)
#        except Exception, e:
#            print _cls(self), ": computing '{}' ({})".format(relpath, e)
#            if callable(f):
#                y = f(*args, **kwargs)
#            else:
#                y = f
#            self._pickle(y, relpath)
#        return y
#
#    def _cached_kernel(self, K, num, relpath, *args, **kwargs):
#        """Loads a kernel from disk, or computes it if necessary.
#
#        :param K: kernel class.
#        :param num: number of entities.
#        :param relpath: path to the Gram matrix relative to the experiment destination directory.
#        :param tol: non-PSD tolerance (default: ``1e-10``).
#        :param do_pairwise: whether to compute the corresponding pairwise kernel (default: ``False).
#        :param pairs: indices of pairs to compute the pairwise kernel for.
#        :param pairwise_op: pairwise operation to performn (default: ``"product"``).
#        """
#        path = os.path.join(self.dst, relpath)
#        force_update = False
#        try:
#            print _cls(self), ": loading '{}'".format(path)
#            assert not self.force_update, "update forced by user"
#            kernel = DummyKernel(path + ".txt", num = num, check_psd = True)
#        except Exception, e:
#            print _cls(self), "|", e
#            print _cls(self), ": computing '{}'".format(path)
#            kernel = K(*args, **kwargs)
#            kernel.check_and_fixup(kwargs.get("tol", 1e-10))
#            kernel.save(path + ".txt")
#            kernel.draw(path + ".png")
#            force_update = True
#
#        do_pairwise = kwargs.get("do_pairwise", False)
#        if not do_pairwise:
#            return kernel
#
#        pairwise_path = path + "-pairwise"
#        try:
#            print _cls(self), ": loading '{}'".format(pairwise_path)
#            assert not self.force_update, "update forced by user"
#            assert not force_update, "sub-kernel updated, forcing update"
#            pairwise_kernel = DummyKernel(pairwise_path + ".txt", check_psd = True)
#        except Exception, e:
#            print _cls(self), "|", e
#            print _cls(self), ": computing '{}'".format(pairwise_path)
#            pairwise_kernel = PairwiseKernel(kwargs.get("pairs"), kernel,
#                                             op = kwargs.get("pairwise_op", "product"))
#            pairwise_kernel.check_and_fixup(kwargs.get("tol", 1e-10))
#            pairwise_kernel.save(pairwise_path + ".txt")
#            pairwise_kernel.draw(pairwise_path + ".png")
#
#        return kernel, pairwise_kernel
#
#    def _load_kernel(self, relpath, num = None, force_update = False):
#        path = os.path.join(self.dst, relpath)
#        print _cls(self), ": loading '{}'".format(path)
#        return DummyKernel(path + ".txt", num = num, check_psd = False)
#
#    def _compute_average_kernel(self, relpaths):
#        matrix = self._load_kernel(relpaths[0]).compute()
#        for relpath in relpaths[1:]:
#            matrix += self._load_kernel(relpath).compute()
#        return matrix * (1.0 / len(relpaths))

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
    def _eval_perf(test_ys, pred_ys):
        fpr, tpr, _ = roc_curve(test_ys, pred_ys)
        pr, rc, f1, sup = precision_recall_fscore_support(test_ys, pred_ys)
        return sup, f1, pr, rc, auc(fpr, tpr)

    def eval_svm(self, folds, ys, kernel, c = 1.0):
        results = []
        for k, (test_indices, train_indices) in enumerate(folds):
            gram = kernel.compute()

            # Split the labels between training and test, and compute the
            # per-class costs on the training labels
            train_ys, test_ys = self._split_vector(ys,
                                                   train_indices,
                                                   test_indices)
            train_gram, test_gram = self._split_matrix(gram,
                                                       train_indices,
                                                       test_indices,
                                                       mode = "sklearn")
            model = SVC(C = c, kernel = "precomputed", class_weight = "auto")
            model.fit(train_gram, train_ys)
            pred_ys = model.predict(test_gram)
            results.append(self._eval_perf(test_ys, pred_ys))
        return results

    def eval_mkl(self, folds, ys, kernels, c = 1.0):
        raise NotImplementedError

    def _compute_kernel(self, Kernel, *args, **kwargs):
        kernel = Kernel(*args, **kwargs)
        kernel.check_and_fixup(kwargs.get("tol", 1e-10))
        return kernel.compute(),

    def _make(self, stages, targets, context={}):
        """Make-like functionality based on "stages"."""

        # Map dependencies to stages
        target_to_stage = {}
        for stage in stages:
            for target in stage.outputs:
                assert not target in target_to_stage, \
                       "two stages produce the same target '{}'".format(target)
                target_to_stage[target] = stage

        def resolve(target_to_stage, target, context, force_update):
            stage = target_to_stage[target]

            # Try to load from the cache (avoids to recur deeper in the
            # dependency graph if all dependencies are satisfied)
            ret, loaded_all = {}, True
            for target in stage.outputs:
                relpath = target + ".pickle"
                try:
                    assert not force_update, "forced update"
                    print _cls(self), ": loading '{}'".format(relpath)
                    ret[target] = self._depickle(relpath)
                except Exception, e:
                    print _cls(self), ": failed to load '{}' ({}), recursing deeper...".format(relpath, e)
                    loaded_all = False
                    break

            if loaded_all:
                return ret
            del ret

            # Resolve the dependencies
            print _cls(self), ": checking for {}".format(stage.inputs)
            for target in stage.inputs:
                if not target in context:

                    # Resolve the current dependency
                    outputs = resolve(target_to_stage, target, context, force_update)

                    # Update the context
                    for target in outputs:
                        assert not target in context
                    context.update(outputs)

            # Now that all dependencies are satisfied, compute the target
            print _cls(self), ": about to run {}".format(stage.f)
            results = stage.f(*[context[target] for target in stage.inputs])
            assert len(results) == len(stage.outputs), \
                   "declared and actual outputs differ: {} vs {}".format(len(results), len(stage.outputs))

            # Prepare the results dictionary and cache them
            ret = {}
            for target, result in zip(stage.outputs, results):
                relpath = target + ".pickle"
                print _cls(self), ": saving '{}'".format(relpath)
                self._pickle(result, relpath)
                ret[target] = result

            return ret

        # Resolve for all targets
        for target in targets:
            resolve(target_to_stage, target, context, self.force_update)
        return context

    def run(self):
        raise NotImplementedError()

from yip09 import *
from sgd import *
from cafa13 import *
from full import *
