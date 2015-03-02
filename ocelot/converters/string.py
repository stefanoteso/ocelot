# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv
from ocelot.services import _cls

class STRINGConverter(Converter):
    """Converter for `STRING <http://string-db.org/>`_.

    .. warning::

        This converter is a STUB.

    :param taxon: list of taxa to process (default: ``"4932"`` aka S. Cerevisiae).
    :param version: STRING DB version to use (default: ``"9.1"``)
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("aliases",         self._siphon_aliases),
            ("interactions",    self._siphon_interactions),
            ("cog",             self._siphon_cog),
        )
        super(STRINGConverter, self).__init__(SUBTARGETS, *args, **kwargs)
        self._taxon = kwargs.get("taxon", "4932")
        self._version = kwargs.get("version", "9.1")

    def _get_path(self, basename):
        import os
        return os.path.join(self.src, "STRING", basename)

    def _siphon_aliases(self, triples):
        FIELDS = ("TAXON", "STRING_ID", "ALIAS_ID", "ALIAS_TYPE")
        path = self._get_path("{}.protein.aliases.v{}.txt") \
            .format(self._taxon, self._version)
        for row in iterate_csv(path, delimiter = "\t", fieldnames = FIELDS,
                               num_skip = 1):
            string_id = O.uri(O.STRING_ID, row["STRING_ID"])
            triples.append((string_id, O.RDF.type, O.OCELOT.STRING_ID))

            alias_id = row["ALIAS_ID"]
            for alias_type in row["ALIAS_TYPE"].split():
                if alias_type == "SGD" and \
                   (alias_id.startswith("S0") or alias_id.startswith("L0")):
                    db_alias_id = O.uri(O.SGD_ID, alias_id)
                else:
                    # XXX handle UniProt ACs here
                    continue
                triples.append((string_id, O.OWL.sameAs, db_alias_id))

    def _siphon_interactions(self, triples):
        pass

    def _siphon_cog(self, triples):
        FIELDS = ("TAXON.STRING_ID", "START", "STOP", "CLUSTER_ID", "ANNOTATION")
        path = self._get_path("COG.mappings.v{}.txt").format(self._version)
        for row in iterate_csv(path, delimiter = "\t", fieldnames = FIELDS,
                               num_skip = 1):
            parts = row["TAXON.STRING_ID"].split(".")
            taxon = parts[0]
            string_id = ".".join(parts[1:])
            if taxon != self._taxon:
                continue
            string_id = O.uri(O.STRING_ID, string_id)
            cluster_id = row["CLUSTER_ID"]
            if cluster_id.startswith("COG"):
                cluster_id = O.uri(O.COG_CLUSTER_ID, cluster_id)
                triples.extend([
                    (string_id, O.STRING_ID_IN_COG, cluster_id),
                    (cluster_id, O.RDF.type, O.COG_CLUSTER_ID)
                ])
            elif cluster_id.startswith("KOG"):
                cluster_id = O.uri(O.KOG_CLUSTER_ID, cluster_id)
                triples.extend([
                    (string_id, O.STRING_ID_IN_KOG, cluster_id),
                    (cluster_id, O.RDF.type, O.COG_CLUSTER_ID)
                ])
            elif cluster_id.startswith("NOG"):
                cluster_id = O.uri(O.NOG_CLUSTER_ID, cluster_id)
                triples.extend([
                    (string_id, O.STRING_ID_IN_NOG, cluster_id),
                    (cluster_id, O.RDF.type, O.COG_CLUSTER_ID)
                ])
