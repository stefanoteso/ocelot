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

    Returns
    -------
    folds : list
        A partition of the proteins as a list of sets of protein IDs.
    """
    ps = permute(ps, rng=rng)
    num_ps_per_fold = len(ps) // num_folds

    folds, base = [], 0
    for k in range(num_folds):
        ps_in_fold = set(ps[base:base + num_ps_per_fold])
        base += num_ps_per_fold
        if k == num_folds - 1:
            ps_in_fold.update(ps[base:])
        folds.append(ps_in_fold)
    assert sum(map(len, folds)) == len(ps)

    return folds

def _distribute_pps(folds, pps, pp_to_feats, state, rng):

    # Map from individual features to pairs
    feat_to_pps = defaultdict(set)
    for pp, feats in pp_to_feats.iteritems():
        if not len(feats):
            feat_to_pps["unannotated"].add(pp)
        else:
            for feat in feats:
                feat_to_pps[feat].add(pp)

    # Distribute the interactions among the folds
    while len(feat_to_pps):

        # Shuffle the folds
        folds = permute(folds, rng)

        # Pick the term with the least unallocated pairs
        cur_feat, symm_pps = min(feat_to_pps.iteritems(),
                                    key=lambda id_pps: len(id_pps[-1]))

        print "distributing interactions with term {} (# interactions = {})".format(cur_feat, len(symm_pps))

        # Desymmetrize the pps to add
        asymm_pps = set()
        for p, q in symm_pps:
            if not (q, p) in asymm_pps:
                asymm_pps.add((p, q))

        # Evenly distribute the associated pairs among all folds.
        for i, (p, q) in enumerate(sorted(asymm_pps)):
            folds[i % len(folds)].update([(p, q, state), (q, p, state)])

        # Update term_to_pps by removing all interactions that we
        # just allocated
        new_feat_to_pps = {}
        for feat in feat_to_pps:
            # Skip the term we just processed (it's done, right?)
            if feat == cur_feat:
                continue
            new_pps = set(pp for pp in feat_to_pps[feat]
                          if not pp in symm_pps)
            if not len(new_pps):
                continue
            new_feat_to_pps[feat] = new_pps
        feat_to_pps = new_feat_to_pps

def _check_ppi_folds(ps, pps):
    assert all((q, p) in pps for (p, q) in pps), \
        "pairs are not symmetric"
    assert all(p in ps for p, _ in pps), \
        "singletons and pairs do not match"

def compute_ppi_folds(ps, pp_pos, pp_neg, p_to_feats, num_folds=10, rng=None):
    """Generates interaction-based folds.

    Folds are composed of *pairs* of (known interacting or assumed
    non-interacting) proteins. Guarantees that:

    - Protein pairs are split evenly according to the feature vectors.
    - The same pair can not occur in distinct folds.
    - The same protein *may* occur in distinct folds.
    - Some proteins may not appear in any fold.

    Parameters
    ----------
    ps : list
        Ordered list of proteins.
    pp_pos : set
        Positive interactions.
    pp_neg : set
        Negative interactions.
    p_to_feats : dict
        Map between proteins to features.
    num_folds : int, optional. (defaults to 10)
        Number of folds.

    Returns
    -------
    folds : list
        Each fold is a set of triples of the form (protein, protein, state).
    """
    folds = [set() for _ in range(num_folds)]

    pp_to_feats = {}
    for p, q in pp_pos | pp_neg:
        p_feats = p_to_feats.get(p, set())
        q_feats = p_to_feats.get(q, set())
        pp_to_feats[(p, q)] = p_feats | q_feats

    _distribute_pps(folds, pp_pos, pp_to_feats, True, rng)
    _distribute_pps(folds, pp_neg, pp_to_feats, False, rng)
    _check_folds(ps, pp_pos | pp_neg, folds)
    return folds
