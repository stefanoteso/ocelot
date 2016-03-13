# -*- coding: utf-8 -*-

import os
import numpy as np
from SPARQLWrapper import SPARQLWrapper, JSON
from sklearn.utils import check_random_state
from sklearn.svm import SVC
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc

import .ontology as O
from . import split_tr_ts

Stage = namedtuple("Stage", ["f", "inputs", "outputs"])

class Scheduler(object):
    """A simple make-like scheduler.

    Parameters
    ----------
    dst : str
        Directory where the results will be put.
    """
    def __init__(self, dst):
        self.dst = dst

    def _resolve_save(self, basename, what):
        """Saves something as basename.pickle in self.dst."""
        relpath = os.path.join(self.dst, basename) + ".pickle"
        print _cls(self), ": saving '{}'".format(relpath)
        try:
            with open(relpath, "wb") as fp:
                pickle.dump(what, fp, protocol=2)
        except Exception, e:
            raise IOError("can not save '{}':\n{}".format(relpath, e))

    def _resolve_load(self, basename):
        """Loads basename.pickle from self.dst."""
        relpath = os.path.join(self.dst, basename) + ".pickle"
        try:
            with open(relpath, "rb") as fp:
                return pickle.load(fp)
        except:
            pass
        raise IOError("can not load '{}'".format(relpath))

    def _resolve(self, target_to_stage, target, context, force_update):
        """Resolves a stage.

        That is, it runs the stage by first making sure that all its
        dependencies are satisfied, recursively.

        Parameters
        ----------
        target_to_stage : dict
            Map from target to Stage.
        target : str
            Target to resolve.
        context : dict
            Map from dependency to value.
        force_update: bool
            Whether to ignore the cache contents.

        Returns
        -------
        ret : dict
            Map from dependency to value, computed by the target stage.
        """
        stage = target_to_stage[target]

        # Try to load the target from the cache. If all dependencies are
        # cached, we do not want to recurse deeper in the dependency graph
        if not force_update:
            ret, all_loaded = {}, True
            for output in stage.outputs:
                try:
                    ret[output] = self._resolve_load(output)
                except Exception, e:
                    print _cls(self), ": failed to load dependency ({}), recursing".format(e)
                    all_loaded = False
            if all_loaded:
                return ret

        # Resolve all dependencies
        print _cls(self), ": resolving deps for '{}'".format(target)
        for input_ in stage.inputs:
            if not input_ in context:
                outputs = self._resolve(target_to_stage, input_, context, force_update)
                for input_ in outputs:
                    assert not input_ in context
                context.update(outputs)

        # Now that all dependencies are satisfied, compute the target
        print _cls(self), ": about to run '{}'".format(stage.f)
        results = stage.f(*[context[input_] for input_ in stage.inputs])
        assert len(results) == len(stage.outputs), \
            "declared and actual outputs differ: {} vs {}".format(len(results), len(stage.outputs))

        # Fill the results dictionary and cache them
        ret = {}
        for output, result in zip(stage.outputs, results):
            self._resolve_save(output, result)
            ret[output] = result

        return ret

    def run(self, stages, targets, context={}, force_update=False):
        """Runs the whole thing."""
        # Map dependencies to stages
        target_to_stage = {}
        for stage in stages:
            for target in stage.outputs:
                assert not target in target_to_stage, \
                    "two stages produce the same target '{}'".format(target)
                target_to_stage[target] = stage

        # Resolve for all targets
        for target in targets:
            self._resolve(target_to_stage, target, context, force_update)

        return context

class Endpoint(object):
    """A simple SPARQL endpoint.

    Parameters
    ----------
    uri : str
        URI of the SPARQL endpoint.
    graph : str
        URI of the default graph.
    """
    def __init__(self, uri, graph):
        self.endpoint, self.graph = SPARQLWrapper(uri), graph
        ans = self.query(u"ASK WHERE {{ GRAPH <{default_graph}> {{ ?s ?p ?o }} }}")
        if not ans[u"boolean"]:
            raise ValueError("no graph '{}' in '{}'".format(graph, uri))

    def query(self, query):
        """Performs a query at the given endpoint.

        Parameters
        ----------
        query : str
            SPARQL query.

        Returns
        -------
        ans : list
            List of triples.
        """
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

class Experiment(object):
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
