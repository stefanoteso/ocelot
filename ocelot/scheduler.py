# -*- coding: utf-8 -*-

import os
import cPickle as pickle
from collections import namedtuple

from ocelot.services import _cls

Stage = namedtuple("Stage", ["f", "inputs", "outputs"])

class Scheduler(object):
    """A simple make-like scheduler."""

    def __init__(self, dst):
        self.dst = dst

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
                return np.load(fp)
        except:
            pass
        raise IOError("can not load '{}'".format(relpath))

    def _resolve(self, target_to_stage, target, context, force_update):
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
