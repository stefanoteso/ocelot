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

    def _resolve_save(self, basename, what):
        relpath = os.path.join(self.dst, basename)
        print _cls(self), ": saving '{}'".format(relpath)
        try:
            what.dump(relpath + ".npy")
            return
        except:
            # Not a numpy object
            pass
        try:
            with open(relpath + ".pickle", "wb") as fp:
                pickle.dump(what, fp)
            return
        except:
            pass
        raise IOError("can not save '{}'".format(relpath))

    def _resolve_load(self, basename):
        relpath = os.path.join(self.dst, basename)
        # Try to load a pickled file
        try:
            with open(relpath + ".pickle", "rb") as fp:
                return pickle.load(fp)
        except:
            pass
        # Try to load a numpy array
        try:
            with open(relpath + ".npy", "rb") as fp:
                return pickle.load(fp)
        except:
            pass
        raise IOError("can not load '{}'".format(relpath))

    def _resolve(self, target_to_stage, target, context, force_update):
        stage = target_to_stage[target]

        # Try to load from the cache: if all dependencies are cached, we
        # do not want to recurse deeper in the dependency graph
        if not force_update:
            ret, loaded_all = {}, True
            for output in stage.outputs:
                try:
                    ret[output] = self._resolve_load(output)
                except Exception, e:
                    print _cls(self), ": failed to load dependency ({}), recursing".format(e)
                    loaded_all = False
            if loaded_all:
                return ret

        # Resolve the dependencies
        print _cls(self), ": resolving dependencies for {}".format(target)
        for input_ in stage.inputs:
            if not input_ in context:

                # Resolve the current dependency
                outputs = self._resolve(target_to_stage, input_, context, force_update)

                # Update the context
                for input_ in outputs:
                    assert not input_ in context
                context.update(outputs)

        # Now that all dependencies are satisfied, compute the target
        print _cls(self), ": about to run {}".format(stage.f)
        results = stage.f(*[context[input_] for input_ in stage.inputs])
        assert len(results) == len(stage.outputs), \
               "declared and actual outputs differ: {} vs {}".format(len(results), len(stage.outputs))

        # Prepare the results dictionary and cache them
        ret = {}
        for output, result in zip(stage.outputs, results):
            self._resolve_save(output, result)
            ret[output] = result

        return ret

    def _make(self, stages, targets, context={}):
        """Make-like functionality based on "stages"."""

        # Map dependencies to stages
        target_to_stage = {}
        for stage in stages:
            for target in stage.outputs:
                assert not target in target_to_stage, \
                       "two stages produce the same target '{}'".format(target)
                target_to_stage[target] = stage

        # Resolve for all targets
        for target in targets:
            self._resolve(target_to_stage, target, context, self.force_update)
        return context

    def run(self):
        raise NotImplementedError()

from yip09 import *
from sgd import *
from cafa13 import *
from full import *
