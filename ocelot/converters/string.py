# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv

class STRINGConverter(Converter):
    """Converter for `STRING <http://string-db.org/>`_.

    The converter assumes that the data directory includes a ``STRING``
    directory with a species-specific STRING dup, laid out like this::

        ${TAXON}.protein.aliases.v${VERSION}.txt            : DONE
        ${TAXON}.protein.links.v${VERSION}.txt              : TODO
        ${TAXON}.protein.links.detailed.v${VERSION}.txt     : TODO
        ${TAXON}.protein.actions.v${VERSION}.txt            : TODO
        ${TAXON}.protein.actions.detailed.v${VERSION}.txt   : TODO
        ${TAXON}.protein.sequences.v${VERSION}.fa           : -
        COG.links.v${VERSION}.txt                           : DONE
        COG.links.detailed.v${VERSION}.txt                  : TODO
        COG.mappings.v${VERSION}.txt                        : TODO

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
        FIELDS = ("ITEM_ID_A", "ITEM_ID_B", "MODE", "ACTION", "A_IS_ACTING",
                  "SCORE", "SOURCES", "TRANSFERRED_SOURCES")
        path = self._get_path("{}.protein.actions.detailed.v{}.txt") \
            .format(self._taxon, self._version)
        for row in iterate_csv(path, delimiter = "\t", fieldnames = FIELDS,
                               num_skip = 1):
            id_a, id_b = row["ITEM_ID_A"], row["ITEM_ID_B"]
            if not (id_a.startswith("{}".format(self._taxon)) and \
                    id_b.startswith("{}".format(self._taxon))):
                continue
            id_a = O.uri(O.STRING_ID, id_a.split(".", 1)[1])
            id_b = O.uri(O.STRING_ID, id_b.split(".", 1)[1])
            mode = O.uri(O.STRING_ACTION_MODE, row["MODE"])
            triples.append((id_a, mode, id_b))

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
