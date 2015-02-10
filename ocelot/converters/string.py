# -*- coding: utf8 -*-

from ocelot.converters import base

class STRINGConverter(base.Converter):

	def __init__(self, *args, **kwargs):
		super(STRINGConverter, self).__init__(*args, **kwargs)
		self.basename = "STRING"
		self.targets = (
			("aliases", self.get_aliases),
			#("interactions", self.get_interactions),
		)

		self.taxa_to_keep = None
		if "taxa_to_keep" in kwargs:
			self.taxa_to_keep = kwargs["taxa_to_keep"]

	def get_path(self, basename):
		return os.path.join(self.src, self.basename, basename)

	def get_aliases(self, triples):

		@unique
		class C(IntEnum):
			TAXON = 0
			STRING_ID = 1
			ALIAS_ID = 2
			ALIAS_TYPE = 3

		for row in self.read_csv(self.get_path("protein.aliases.v9.1.taxon=4932.txt"), C, skip = "#"):

			if not "UniProt_AC" in row[C.ALIAS_TYPE]:
				continue

			if self.taxa_to_keep != None and not row[C.TAXON] in self.taxa_to_keep:
				continue

			string_id	= OSBR.STRING_ID_ROOT.value + sanitize(row[C.STRING_ID])
			alias_id	= OSBR.UNIPROT_ID_ROOT.value + sanitize(row[C.ALIAS_ID])

			triples.extend([
				(string_id, NS.RDF.type, OSBR.STRING_ID.value),
				(string_id, NS.OWL.sameAs, alias_id),
			])
