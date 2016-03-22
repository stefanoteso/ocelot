# -*- coding: utf-8 -*-

import os
import cPickle as pickle
import numpy as np
from collections import namedtuple, defaultdict
from SPARQLWrapper import SPARQLWrapper, JSON
from sklearn.utils import check_random_state

from ocelot import *
from ocelot.ontology import BINDINGS
from ocelot.utils import permute

Stage = namedtuple("Stage", ["f", "inputs", "outputs"])

PHONY = lambda *args, **kwargs: (None,)

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
        print "saving '{}'".format(relpath)
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
                    print "failed to load dependency ({}), recursing".format(e)
                    all_loaded = False
            if all_loaded:
                return ret

        # Resolve all dependencies
        print "resolving deps for '{}'".format(target)
        for input_ in stage.inputs:
            if not input_ in context:
                outputs = self._resolve(input_, context, force)
                for input_ in outputs:
                    assert not input_ in context
                context.update(outputs)

        # Now that all dependencies are satisfied, compute the target
        print "about to run '{}'".format(stage)
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
        for shortcut, namespace in BINDINGS:
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

    def compute_kernel(self, Kernel, png_path, *args, **kwargs):
        """Wrapper to compute a kernel and fix it up."""
        # TODO move the logic to ocelot.Kernel
        tol = kwargs.pop("tol", 1e-10)
        kernel = Kernel(*args, **kwargs)
        kernel.check_and_fixup(tol)
        kernel.draw(os.path.join(self.dst, png_path))
        return kernel.compute(),


def compute_p_folds(ps, p_to_features, num_folds=10, rng=None):
    """Generates balanced protein folds.

    Parameters
    ----------
    ps : list
        Ordered collection of protein IDs.
    p_to_features : dict
        Map from protein IDs to sets of features.
    num_folds : int, optional. (defaults to 10)
        Number of folds to generate.
    rng : None or int or numpy.random.RandomStream, optional. (defaults to None)
        The RNG.

    Returns
    -------
    folds : list
        A partition of the proteins as a list of sets of protein IDs.
    """
    # Map features to proteins that have them
    feature_to_ps = defaultdict(set)
    for p, features in p_to_features.iteritems():
        if not len(features):
            feature_to_ps["unannotated"].add(p)
        else:
            for feature in features:
                feature_to_ps[feature].add(p)

    # Build the folds
    folds = [set() for _ in range(num_folds)]
    while len(feature_to_ps):

        # Pick the feature with the least unallocated proteins
        best_feature, best_ps = min(feature_to_ps.iteritems(),
                                    key=lambda feature_ps: len(feature_ps[-1]))

        print "picked feature", best_feature, ", distributing", len(best_ps), "proteins"

        # Evenly distribute the associated proteins among all folds
        folds = permute(folds, rng)
        for i, p in enumerate(sorted(best_ps)):
            folds[i % num_folds].add(p)

        # Update the map
        new_feature_to_ps = {}
        for feature, ps in feature_to_ps.iteritems():
            if feature == best_feature:
                continue
            unallocated_ps = set(ps) - set(best_ps)
            if len(unallocated_ps):
                new_feature_to_ps[feature] = unallocated_ps
        feature_to_ps = new_feature_to_ps

    # Check that the folds make sense
    for p in ps:
        assert sum(int(p in fold) for fold in folds) == 1

    # Compute per-term fold unbalance
    all_features = set()
    for p, features in p_to_features.iteritems():
        all_features.update(features)

    p_to_i = {p: i for i, p in enumerate(ps)}
    feature_to_j = {feature: i for i, feature in enumerate(all_features)}

    phi = np.zeros((num_folds, len(feature_to_j)))
    for k, fold in enumerate(folds):
        for p in fold:
            phi[k, [feature_to_j[feature] for feature in p_to_features[p]]] += 1
    avg_phi = np.mean(phi, axis=0)
    unbalance = np.sum(np.abs(phi - avg_phi), axis=1) / len(feature_to_j)

    return folds, unbalance

def _distribute_pps(p_folds, pps, state, rng):
    pp_folds = defaultdict(set)

    # De-symmetrize the interactions
    asymm_pps = set()
    for p, q in pps:
        if not (q, p) in asymm_pps:
            asymm_pps.add((p, q))

    # Shuffle the interactions
    asymm_pps = permute(sorted(asymm_pps), rng=rng)

    # Assign them to folds sequentially
    base, n = 0, len(asymm_pps) / len(p_folds)
    for k in range(len(p_folds)):
        if k == len(p_folds) - 1:
            pp_folds[k].update({(p, q, state) for p, q in asymm_pps[n*k:]})
        else:
            pp_folds[k].update({(p, q, state) for p, q in asymm_pps[n*k:n*(k + 1)]})

    # Symmetrize the folds
    for pp_fold in pp_folds.itervalues():
        pp_fold.update({(q, p, state) for p, q, state in pp_fold})

    return pp_folds

def compute_pp_folds(p_folds, pp_pos, pp_neg, rng=None):
    pp_folds_pos = _distribute_pps(p_folds, pp_pos, True,  rng)
    pp_folds_neg = _distribute_pps(p_folds, pp_neg, False, rng)

    pp_folds = []
    for k in range(len(p_folds)):
        pp_folds.append(pp_folds_pos[k] | pp_folds_neg[k])

    print map(len, pp_folds)

    for k in range(len(p_folds)):
        num_pps_for_p = 0.0
        for p in p_folds[k]:
            for s, t, state in pp_folds[k]:
                if s == p or t == p:
                    num_pps_for_p += 1
        avg_pps_for_p = num_pps_for_p / len(p_folds[k])
        print avg_pps_for_p

    for pp_fold in pp_folds:
        assert all((q, p, state) in pp_fold for p, q, state in pp_fold), \
            "pp folds are not symmetric"

    return pp_folds
