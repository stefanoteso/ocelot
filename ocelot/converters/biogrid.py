# -*- coding: utf8 -*-

from ocelot.converters import base

class BioGRIDConverter(base.Converter):
	"""
	BioGRID PSI-MI Tab 2.5 converter.
	"""

	def __init__(self, *args, **kwargs):
		super(BioGRIDConverter, self).__init__(*args, **kwargs)
		self.path = kwargs["path"]
		self.basename = "biogrid"
		self.targets = (
			("interactions", self.get_interactions),
			("xrefs", self.get_xrefs),
		)

	def get_xrefs(self, triples):

		@unique
		class C(IntEnum):
			# Identifier for this ID within the BioGRID Database
			BIOGRID_ID			= 0
			# Identifier Value Mapped
			IDENTIFIER_VALUE		= 1
			# A Brief Text Description of the Identifier
			IDENTIFIER_TYPE			= 2
			# Official name of the organism
			ORGANISM_OFFICIAL_NAME		= 3

		print "XXX get_xrefs()"

	def get_interactions(self, triples):

		@unique
		class C(IntEnum):
			# Unique identifier for interactor [A,B]. Represented as
			# databaseName:id. Here databaseName is defined by the PSI-MI
			# vocabulary. Can be multiple.
			ID_INTERACTOR_A			= 0
			ID_INTERACTOR_B			= 1
			# Alt. identifiers for interactor [AB]. Same format as above.
			ALT_IDS_INTERACTOR_A		= 2
			ALT_IDS_INTERACTOR_B		= 3
			# Aliases for interactor [AB]. Same format as above.
			ALIASES_INTERACTOR_A		= 4
			ALIASES_INTERACTOR_B		= 5
			# Interaction detection method taken from the PSI-MI vocabulary.
			# Represented as "databaseName:identifier(methodName)". Can be
			# multiple.
			INTERACTION_DETECTION_METHOD	= 6
			# What it says. Can be multiple.
			PUBLICATION_1ST_AUTHOR		= 7
			# What it says. Can be multiple.
			PUBLICATION_IDENTIFIERS		= 8
			# NCBI Taxonomy ID for interactor [AB]. The database name for NCBI
			# taxid is taken from the PSI-MI vocabulary. Represented as
			# "databaseName:identifier" or "databaseName:identifier(species)"
			# (typically databaseName is "taxid"). Can be multiple.
			#
			# Special values are "-1" (in vitro), "-2" (chemical synthesis),
			# "-3" (unknown), "-4" (in vivo), "-5" (in silico).
			TAXID_INTERACTOR_A		= 9
			TAXID_INTERACTOR_B		= 10
			# Interaction types taken from the PSI-MI vocabulary. Represented
			# as "databaseName:identifier(interactionType)". Can be multiple.
			INTERACTION_TYPES		= 11
			# Source database taken from the PSI-MI vocabulary. Represented
			# as "databaseName:identifier(sourceName)" Can be multiple.
			SOURCE_DATABASE			= 12
			# Interaction identifier(s) in the corresponding source database.
			# Represented as "databaseName:identifier".
			INTERACTION_IDENTIFIERS		= 13
			# Confidence score that is *not* taken from PSI-MI vocabulary. Can
			# be multiple.
			CONFIDENCE_VALUES		= 14

		def parse_id_list(s):
			"""
			Given a string like

				'foo:a|bar:c|foo:b'
				
			returns a dictionary like:

				{ 'foo' : ['a','b'], 'bar' : ['c'] }
			"""
			if s == '-':
				return None
			results = {}
			for word in s.split("|"):
				# there are IDs like:     entrez gene/locuslink:"zgc:136271"
				parts = word.split(":")
				if len(parts) == 2:
					db, id = parts[0], parts[1]
				elif len(parts) == 3:
					db, id = parts[0], parts[1] + ":" + parts[2]
				else:
					print word
					assert(0)
				if not db in results:
					results[db] = []
				results[db].append(id.replace('"', ''))
				
			return results

		def parse_tax_id(s):
			if s == '-' or not s.startswith("taxid:") or '|' in s:
				return None
			return rdflib.URIRef(NS.TAX + s.split(":")[1])

		def parse_int_id(s):
			if s == '-' or not s.startswith("BIOGRID:") or '|' in s:
				return None
			return rdflib.URIRef(NS.BIOGRID + "interaction_" + s.split(":")[1])

		def parse_psimi_terms(s):
			if s == '-' or not "psi-mi" in s:
				return None
			terms = []
			for part in s.split("|"):
				#
				# the text looks like this:
				#
				#     psi-mi:"MI:0915"(physical association)
				#
				lo_dquotes = s.find('"')
				hi_dquotes = s.rfind('"')
				assert(hi_dquotes - lo_dquotes == 8)
				term = part[(lo_dquotes + 1) : hi_dquotes].replace(":", "_")
				terms.append(rdflib.URIRef(NS.OBO + term))
			return terms

		for row in self.read_csv(os.path.join(self.src, self.path), C, skip = '#'):

			ids_a     = parse_id_list(row[C.ID_INTERACTOR_A])
			alt_ids_a = parse_id_list(row[C.ALT_IDS_INTERACTOR_A])
			aliases_a = parse_id_list(row[C.ALIASES_INTERACTOR_A])

			assert(ids_a != None and len(ids_a["BIOGRID"]) == 1)
			biogrid_id_a = rdflib.URIRef(NS.BIOGRID + ids_a["BIOGRID"][0])

			ids_b     = parse_id_list(row[C.ID_INTERACTOR_B])
			alt_ids_b = parse_id_list(row[C.ALT_IDS_INTERACTOR_B])
			aliases_b = parse_id_list(row[C.ALIASES_INTERACTOR_B])

			assert(ids_b != None and len(ids_b["BIOGRID"]) == 1)
			biogrid_id_b = rdflib.URIRef(NS.BIOGRID + ids_b["BIOGRID"][0])

			taxid_a		= parse_tax_id(row[C.TAXID_INTERACTOR_A])
			taxid_b		= parse_tax_id(row[C.TAXID_INTERACTOR_B])

			interaction	= parse_int_id(row[C.INTERACTION_IDENTIFIERS])
			int_methods	= parse_psimi_terms(row[C.INTERACTION_DETECTION_METHOD])
			int_types	= parse_psimi_terms(row[C.INTERACTION_TYPES])
			int_sources	= parse_psimi_terms(row[C.SOURCE_DATABASE])

			assert(len(int_methods) == 1)
			assert(len(int_types) == 1)
			assert(len(int_sources) == 1)

			triples.extend([
				(interaction, NS.RDF.type, NS.BIOGRID.interaction),
				(interaction, NS.BIOGRID.interactor, biogrid_id_a),
				(interaction, NS.BIOGRID.interactor, biogrid_id_b),
				(biogrid_id_a, NS.BIOGRID.taxid, taxid_a),
				(biogrid_id_b, NS.BIOGRID.taxid, taxid_b),
			])

			for uri in int_methods:
				triples.append((interaction, NS.BIOGRID.int_method, uri))
			for uri in int_types:
				triples.append((interaction, NS.BIOGRID.int_type, uri))
			for uri in int_sources:
				triples.append((interaction, NS.BIOGRID.int_source, uri))
