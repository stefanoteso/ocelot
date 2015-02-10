# -*- coding: utf8 -*-

import SPARQLWrapper as sw
import networkx as nx

PREFIXES = """
PREFIX osbr: <http://ontosbr.org/>
PREFIX pdbo: <http://rdf.wwpdb.org/schema/pdbx-v40.owl#> 
"""

def _valueof(bindings, varname, default = None):
	if not varname in bindings:
		return default
	# TODO add support for literals
	uri = bindings[varname]["value"]
	return uri.split("/")[-1]

class Endpoint(object):
	"""A wrapper around a SPARQLWrapper endpoint."""

	def __init__(self, endpoint, debug = False):
		self.wrapper = sw.SPARQLWrapper(endpoint)
		self.debug = debug

	def _print_bindings(self, bindings):
		for binding in bindings:
			pprint(sorted([(key, binding[key]["value"]) for key in binding]))

	def query(self, query):
		query = "{}\n{}".format(PREFIXES, query)
		if self.debug:
			print "DEBUG: executing query:\b{}".format(query)
		self.wrapper.setQuery(query)
		self.wrapper.setReturnFormat(sw.JSON)
		bindings = self.wrapper.query().convert()["results"]["bindings"]
		if self.debug:
			print "DEBUG: query returned:"
			self._print_bindings(bindings)
		return bindings

class PINBuilder(object):
	"""
	Class to build a networkx PIN out of a OSBR virtuoso graph.

	XXX we currently treat interactions as undirected, while SGD provides
	information about which feature is the bait and which is the hit.
	"""

	def __init__(self, graph, endpoint, debug = False):
		self.endpoint = Endpoint(endpoint, debug = debug)
		self.graph = graph
		self.debug = debug

	def build(self):

		pin = nx.Graph()

		print "PIN: fetching proteins"

		query = """
		SELECT DISTINCT ?orf ?seq ?fun
		FROM <""" + self.graph + """>
		WHERE {
			?orf a osbr:sgd_id ;
				osbr:sgd_id_has_type osbr:sgd_feature_type.ORF ;
				osbr:sgd_id_has_qualifier osbr:sgd_id_qualifier.Verified .

			OPTIONAL { ?orf osbr:sgd_id_has_sequence ?seq . }
			OPTIONAL { ?orf osbr:sgd_id_has_goslim ?fun . }
		}
		"""

		num = 0
		for bindings in self.endpoint.query(query):
			orf = _valueof(bindings, "orf")
			seq = _valueof(bindings, "seq")
			fun = _valueof(bindings, "fun")

			assert(orf)
			assert(seq)

			pin.add_node(orf, seq = seq)
			if fun:
				if not "fun" in pin.node[orf]:
					pin.node[orf]["fun"] = set()
				pin.node[orf]["fun"].add(fun)

			num += 1

		print "PIN: found {} proteins".format(num)

		print "PIN: fetching protein-protein interactions"

		query = """
		SELECT DISTINCT ?b ?h
		WHERE {
			?b a osbr:sgd_id ;
				owl:sameAs ?b_feature ;
				osbr:sgd_id_has_type osbr:sgd_feature_type.ORF ;
				osbr:sgd_id_has_qualifier osbr:sgd_id_qualifier.Verified .

			?h a osbr:sgd_id ;
				owl:sameAs ?h_feature ;
				osbr:sgd_id_has_type osbr:sgd_feature_type.ORF ;
				osbr:sgd_id_has_qualifier osbr:sgd_id_qualifier.Verified .

			?interaction a osbr:sgd_int ;
				osbr:sgd_int_has_bait ?b_feature ;
				osbr:sgd_int_has_hit ?h_feature ;
				osbr:sgd_int_has_type osbr:sgd_int_type.physical_interactions ;
				osbr:sgd_int_has_cur_type osbr:sgd_int_cur_type.manually_curated .
		}
		"""

		num = 0
		for bindings in self.endpoint.query(query):
			b	= _valueof(bindings, "b")
			h	= _valueof(bindings, "h")

			pin.add_edge(b, h)

			num += 1

		print "PIN: found {} interactions".format(num)

		return pin

class DINBuilder(object):
	"""Converts Ocelot RDF data into a domain-domain interaction network."""

	def __init__(self, graph, endpoint, debug = False):
		self.endpoint = Endpoint(endpoint, debug = debug)
		self.graph = graph
		self.debug = debug

	def build(self):
		din = nx.Graph()

		print "DIN: fetching domain instances"

		query = """
		SELECT DISTINCT ?region ?pfam ?struct
		FROM <""" + self.graph + """>
		WHERE {
			?region a osbr:ipfam_region ;
				osbr:ipfam_region_instance_of ?pfam ;
				osbr:ipfam_region_occurs_in ?struct .

			?pfam a osbr:pfam_id ;
				osbr:pfam_id_has_type osbr:pfam_type.family .
		}
		"""

		num = 0
		for bindings in self.endpoint.query(query):
			region	= _valueof(bindings, "region")
			pfam	= _valueof(bindings, "pfam")
			struct	= _valueof(bindings, "struct")

			#print region, pfam, struct (pdb_id + chain)

			num += 1

		print "DIN: found {} domain instances".format(num)

		print "DIN: fetching domain instance-instance interactions"

		query = """
		SELECT DISTINCT ?region_int ?complex
		FROM <""" + self.graph + """>
		WHERE {
			?region_int a osbr:ipfam_region_int ;
				osbr:ipfam_region_int_occurs_in ?complex .

			?region a osbr:ipfam_region ;
				osbr:ipfam_region_in_region_int ?region_int .
		}
		"""

		num = 0
		for bindings in self.endpoint.query(query):
			region_int	= _valueof(bindings, "region_int")
			complex		= _valueof(bindings, "complex")

			num += 1

		print "DIN: found {} domain instance-instance interactions".format(num)

		return din

class RINBuilder(object):
	"""Converts Ocelot RDF data into a residue-residue interaction network."""

	def __init__(self, graph, endpoint, debug = False):
		self.endpoint = Endpoint(endpoint, debug = debug)
		self.graph = graph
		self.debug = debug

	def build(self):
		rin = nx.Graph()
		return rin

class NegPINBuilder(object):
	pass

class ProteinKernelBuilder(object):
	pass

class DomainKernelBuilder(object):
	pass

class ResidueKernelBuilder(object):
	pass

