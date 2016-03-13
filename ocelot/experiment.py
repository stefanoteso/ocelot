# -*- coding: utf-8 -*-

import os
from collections import namedtuple
from SPARQLWrapper import SPARQLWrapper, JSON
from sklearn.utils import check_random_state
from sklearn.svm import SVC
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc

import ontology as O
from ocelot import *
from ocelot.services import _cls

Stage = namedtuple("Stage", ["f", "inputs", "outputs"])

PHONY = lambda *args, **kwargs: None,

class Scheduler(object):
    """A simple make-like scheduler.

    Parameters
    ----------
    stages : list
        List of Stage objects.
    dst : str
        Directory where the results will be put.
    """
    def __init__(self, stages, dst):
        target_to_stage = {}
        for stage in stages:
            assert len(stage.outputs), \
                "stage '{}' produces no output".format(stage)
            for target in stage.outputs:
                assert not target in target_to_stage, \
                    "two stages produce the same target '{}'".format(target)
                target_to_stage[target] = stage
        self._target_to_stage = target_to_stage
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

    def _resolve(self, target, context, force):
        """Resolves a stage.

        It runs each target by first making sure that all its dependencies are
        satisfied, recursively. Dependencies and results are gathered in a
        single big dict, the context.

        Parameters
        ----------
        target : str
            Target to resolve.
        context : dict
            Map from dependency to value.
        force: bool
            Whether to ignore the cache contents.

        Returns
        -------
        ret : dict
            Map from dependency to value, computed by the target stage.
        """
        stage = self._target_to_stage[target]

        # Try to load the target from the cache. If all dependencies are
        # cached, we do not want to recurse deeper in the dependency graph
        if not force:
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
                outputs = self._resolve(input_, context, force)
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

    def run(self, targets, context={}, force=False):
        """Runs the whole thing.

        Parameters
        ----------
        targets : list of str
            List of targets to be resolved.
        context : dict, optional. (defaults to {})
            Map from dependencies to values.
        force : bool, optional. (defaults to False)
            Whether to ignore the cache.

        Returns
        -------
        context : dict
            The context, updated with the values of the targets.
        """
        for target in targets:
            self._resolve(target, context, force)
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
        ans = self.query(u"ASK WHERE {{ GRAPH <{graph}> {{ ?s ?p ?o }} }}")
        if not ans[u"boolean"]:
            raise ValueError("no graph '{}' in '{}'".format(graph, uri))

    @staticmethod
    def _cast(d):
        """Converts SPARQLWrapper values to standard Python values."""
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

    def query(self, query):
        """Performs a query.

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
        query = prefixes + query.format(graph=self.graph)
        self.endpoint.setQuery(query)
        self.endpoint.setReturnFormat(JSON)
        return self.endpoint.query().convert()

    def iterquery(self, query, n=None):
        """Performs a query and yields the parsed results.

        Parameters
        ----------
        query : str
            SPARQL query.

        Returns
        -------
        bindings : dict
            Map from variable name to value.
        """
        ans = self.query(query)
        assert ans and len(ans) and "results" in ans
        for bindings in ans["results"]["bindings"]:
            bindings = {k: self._cast(v) for k, v in bindings.iteritems()}
            if not n is None:
                assert len(bindings) == n, bindings
            yield bindings

class Experiment(object):
    """Base class for all experiments.

    Parameters
    ----------
    stages : list
        List of Stage objects.
    src : str
        Directory where the raw database data is held.
    dst : str
        Directory where the results will be held.
    endpoint : Endpoint
        The SPARQL endpoint.
    rng : np.random.RandomStream or int or None, optional. (defaults to None)
        RNG.
    """
    def __init__(self, stages, src, dst, endpoint, rng=None):
        try:
            os.mkdir(dst)
        except:
            pass
        if not (src and os.path.isdir(src)):
            raise ValueError("'{}' is not a valid directory".format(src))
        if not (dst and os.path.isdir(dst)):
            raise ValueError("'{}' is not a valid directory".format(dst))
        self.src, self.dst = src, dst

        self.endpoint = endpoint
        self._scheduler = Scheduler(stages, self.dst)
        self._rng = check_random_state(rng)

    def run(self, targets=None, context={}, force=False):
        """Executes the targets with the given context.

        Parameters
        ----------
        targets : list, optional. (defaults to ["__all"])
            List of targets to execute.
        context : dict, optional. (defaults to {})
            Context.
        force: bool, optional. (defaults to False)
            Whether to ignore cached results altogether.
        """
        targets = ["__all"] if targets is None else targets
        self._scheduler.run(targets, context=context, force=force)

    @staticmethod
    def _compute_kernel(self, Kernel, *args, **kwargs):
        """Wrapper to compute a kernel and fix it up."""
        tol = kwargs.pop("tol", 1e-10)
        kernel = Kernel(*args, **kwargs)
        kernel.check_and_fixup(tol)
        return kernel.compute(),

    @staticmethod
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
