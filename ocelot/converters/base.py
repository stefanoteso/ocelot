# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.services import iterate_csv
from rdflib import Namespace as NS, URIRef as U, BNode as B, Literal as L, Graph

import os

class Converter(object):
    """Abstract class for database converters.

    :param src: path to the source databases.
    :param dst: path to the destination folder.
    :param basename: name of this target.
    :param path: base path of this target data within ``src``.
    :param force_update: convert data even if the target RDF file exists.
    """

    def __init__(self, all_subtargets, src, dst, **kwargs):
        self.src = src
        self.dst = dst
        self.basename = kwargs["basename"]
        self.path = kwargs.get("path", None)
        self.force_update = kwargs.get("force_update", False)

        known_subtarget_names = [ pair[0] for pair in all_subtargets ]
        requested_subtargets = kwargs.get("subtargets", all_subtargets)

        subtargets = []
        for subtarget, method in requested_subtargets:
            if not subtarget in known_subtarget_names:
                print "ERROR: invalid subtarget '{}:{}', ignoring".format(self.basename, subtarget)
            else:
                subtargets.append((subtarget, method))
        self.subtargets = subtargets

    def siphon(self):
        for name, method in self.subtargets:
            rdf_path = os.path.join(self.dst, "{}-{}.ttl".format(self.basename, name))
            if os.path.isfile(rdf_path) and not self.force_update:
                print "Target '{}' exists, skipping".format(rdf_path)
                continue

            print "Converting {}:{} into '{}'".format(self.basename, name, rdf_path)
            graph = Graph()
            for shortcut, namespace in O.BINDINGS:
                graph.bind(shortcut, namespace)
            triples = []
            method(triples)
            for triple in triples:
                graph.add(triple)
            graph.serialize(destination = rdf_path, format = "turtle")
            del triples
            del graph

    @staticmethod
    def _add_if_new(triples, a_set, s, p, o):
        if not s in a_set:
            a_set.add(s)
            triples.append((s, p, o))

class CSVConverter(Converter):
    """Base class to convert a CSV file into RDF triples.

    :param fields: list of CSV field names.
    :param delimiter: CSV delimiter.
    """
    def __init__(self, fields, delimiter, num_skip = 0, *args, **kwargs):
        SUBTARGETS = (("rows", self._siphon_rows),)
        super(CSVConverter, self).__init__(SUBTARGETS,
                                           *args, **kwargs)
        self.fields = fields
        self.delimiter = delimiter
        self.num_skip = num_skip
    def _siphon_rows(self, triples):
        import csv
        path = os.path.join(self.src, self.path)
        for row in iterate_csv(path, delimiter = self.delimiter,
                               fieldnames = self.fields,
                               num_skip = self.num_skip):
            self._siphon_row(triples, row)
    def _siphon_row(self, triples, row):
        raise NotImplementedError

