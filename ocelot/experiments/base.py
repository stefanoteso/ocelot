# -*- coding: utf-8 -*-

import ocelot.ontology as O

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

    def _crossvalidate_svm(self, ys, k, folds = None):
        """NOTE: y must be {-1,+1}. """
        from modshogun import RealFeatures, BinaryLabels, CustomKernel, LibSVM

        kernel = CustomKernel()
        kernel.set_full_kernel_matrix_from_full(k._matrix)

        labels = BinaryLabels(ys)
        svm = LibSVM(1, kernel, labels)
        svm.train()
        pred_ys = svm.apply().get_labels()

    def _crossvalidate_mkl(self, ys, ks, folds = None):
        from modshogun import CombinedKernel, CustomKernel
        from modshogun import MKLClassification

        combined_kernel = CombinedKernel()
        for k in ks:
            combined_kernel.append_kernel(CustomKernel(ks.compute()))

        model = MKLClassification()
        model.set_mkl_norm(1) # 2, 3
        model.set_C(1, 1) # positive, negative cost
        model.set_kernel(combined_kernel)
        model.set_labels(BinaryLabels(ys))

        model.train()
        subkernel_weights = combined_kernel.get_subkernel_weights()
        print subkernel_weights

        predictions = mkl.apply()
        print predictions

    def _crossvalidate_sbr(self, ys, ks, folds = None):
        pass

    def run(self):
        raise NotImplementedError

