# -*- coding: utf8 -*-

from rdflib import Namespace as NS, URIRef as U, BNode as B, Literal as L, Graph
from pprint import pprint

import ocelot.ontology as O
from ocelot.services import iterate_csv

import os

class Converter(object):
    """Abstract class for database converters."""

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

            print "Converting {}:{} into '{}'".format(self.basename, name, rdf_path)

            graph = Graph()

            for shortcut, namespace in O.BINDINGS:
                graph.bind(shortcut, namespace)

            if os.path.isfile(rdf_path) and not self.force_update:
                print "RDF file '{}' exists, skipping".format(rdf_path)
                continue

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

    :param basename: source basename, e.g. 'mint' or 'intact'.
    :param fields: list of CSV field names.
    """
    def __init__(self, fields, delimiter, *args, **kwargs):
        SUBTARGETS = (("rows", self._siphon_rows),)
        super(CSVConverter, self).__init__(SUBTARGETS,
                                           *args, **kwargs)
        self.fields = fields
        self.delimiter = delimiter
    def _siphon_rows(self, triples):
        import csv
        path = os.path.join(self.src, self.path)
        for row in iterate_csv(path, delimiter = self.delimiter,
                               fieldnames = self.fields):
            self._siphon_row(triples, row)
    def _siphon_row(self, triples, row):
        raise NotImplementedError

