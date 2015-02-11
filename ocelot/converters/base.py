# -*- coding: utf8 -*-

from rdflib import Namespace as NS, URIRef as U, BNode as B, Literal as L, Graph
from pprint import pprint

import ocelot.ontology as O

import os

class Converter(object):
    """Abstract class for database converters."""

    def __init__(self, basename, all_subtargets, src, dst, **kwargs):
        self.basename = basename
        self.src = src
        self.dst = dst
        self.path = kwargs.get("path", None)
        self.force_update = kwargs.get("force_update", False)

        known_subtarget_names = [ pair[0] for pair in all_subtargets ]
        requested_subtargets = kwargs.get("subtargets", all_subtargets)

        subtargets = []
        for subtarget, method in requested_subtargets:
            if not subtarget in known_subtarget_names:
                print "ERROR: invalid subtarget '{}:{}', ignoring".format(basename, subtarget)
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

    def sanitize(self, s):
        REPS = (
            (" ", "_"),
            ("'", "_SQUOTE_"),
            ('"', "_DQUOTE_"),
            ("(", "_"),
            (")", "_")
        )
        for f, t in REPS:
            s = s.replace(f, t)
        return s

    @staticmethod
    def _add_if_new(triples, a_set, s, p, o):
        if not s in a_set:
            a_set.add(s)
            triples.append((s, p, o))

    def pdb_id_to_uri(self, id, chain = None):
        base_uri = id.lower()
        if chain:
            base_uri += "_" + chain.upper()
        return U(NS.PDBR + base_uri)

    def uniprot_id_to_uri(self, id):
        return U(O.UNIPROT_ID_ROOT.value + id())

def iterate_csv(path, skip_first = False, **kwargs):
    import csv
    with open(path, "rU") as fp:
        reader = csv.DictReader(fp, **kwargs)
        for row_as_dict in reader:
            if skip_first:
                skip_first = False
                continue
            yield row_as_dict

class CSVConverter(Converter):
    """Base class to convert a CSV file into RDF triples.

    :param basename: source basename, e.g. 'mint' or 'intact'.
    :param fields: list of CSV field names.
    """
    def __init__(self, basename, fields, delimiter, *args, **kwargs):
        SUBTARGETS = (("rows", self._siphon_rows),)
        super(CSVConverter, self).__init__(basename, SUBTARGETS,
                                           *args, **kwargs)
        self.fields = fields
        self.delimiter = delimiter
    def _siphon_rows(self, triples):
        import csv
        path = os.path.join(self.src, self.path)
        with open(path, "rt") as fp:
            reader = csv.DictReader(fp, fieldnames = self.fields,
                                    restkey = "UNKNOWN",
                                    delimiter = self.delimiter)
            for row in reader:
                self._siphon_row(triples, row)
    def _siphon_row(self, triples, row):
        raise NotImplementedError

