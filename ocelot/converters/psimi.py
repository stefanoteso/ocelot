# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import CSVConverter
from rdflib import URIRef as U, BNode as B, Literal as L

# Columns for PSI-MI TAB versions 2.5, 2.6 and 2.7 [1,2,3]. All columns are
# mandatory.
_FIELDS = (
    # 2.5
    "ID_INTERACTOR_A",
    "ID_INTERACTOR_B",
    "ALT_IDS_INTERACTOR_A",
    "ALT_IDS_INTERACTOR_B",
    "ALIASES_INTERACTOR_A",
    "ALIASES_INTERACTOR_B",
    "INTERACTION_DETECTION_METHOD",
    "PUBLICATION_1ST_AUTHOR",
    "PUBLICATION_IDENTIFIERS",
    "TAXID_INTERACTOR_A",
    "TAXID_INTERACTOR_B",
    "INTERACTION_TYPES",
    "SOURCE_DATABASE",
    "INTERACTION_IDENTIFIERS",
    "CONFIDENCE_VALUES",
    # 2.6
    "COMPLEX_EXPANSION",
    "BIOLOGICAL_ROLE_A",
    "BIOLOGICAL_ROLE_B",
    "EXPERIMENTAL_ROLE_A",
    "EXPERIMENTAL_ROLE_B",
    "INTERACTOR_TYPE_A",
    "INTERACTOR_TYPE_B",
    "XREF_A",
    "XREF_B",
    "XREF_INTERACTION",
    "ANNOTATION_A",
    "ANNOTATION_B",
    "ANNOTATION_INTERACTION",
    "TAXID_HOST",
    "INTERACTION_PARAMETERS",
    "CURATION_CREATION_TIME",
    "CURATION_UPDATE_TIME",
    "CHECKSUM_A",
    "CHECKSUM_B",
    "CHECKSUM_INTERACTION",
    "IS_NEGATIVE",
    # 2.7
    "FEATURES_A",
    "FEATURES_B",
    "STOICHIOMETRY_A",
    "STOICHIOMETRY_B",
    "PARTICIPANT_IDENTIFICATION_METHOD_A",
    "PARTICIPANT_IDENTIFICATION_METHOD_B",
)

class PsiMiTabConverter(CSVConverter):
    """PSI-MI TAB Converter

    Supports the 2.5, 2.6 and 2.7 versions.

    :param basename: source basename, e.g. 'mint' or 'intact'.

    *References*

    .. [mitab25] https://code.google.com/p/psicquic/wiki/MITAB25Format
    .. [mitab26] https://code.google.com/p/psicquic/wiki/MITAB26Format
    .. [mitab27] https://code.google.com/p/psicquic/wiki/MITAB27Format
    """
    def __init__(self, *args, **kwargs):
        super(PsiMiTabConverter, self).__init__(_FIELDS, "\t",
                                                *args, **kwargs)

    @staticmethod
    def _get_parts(string):
        """Converts a pipe-separated PSI-MI string to a list of strings.

        :param string: the string to be split.
        """
        if not string or string == "-":
            return []
        return filter(lambda part: not "unknown" in part, string.split("|"))

    @staticmethod
    def _prefix_to_uri(string):
        """Converts a prefixed string to a URI.

        :param: the string to be de-prefixed.
        """
        PREFIXES = (
            ("psi-mi", O.OBO),
            ("imex", O.OBO.IMEX),
            ("intact", O.OBO.INTACT),
            ("mint", O.OBO.MINT),
            ("uniprotkb", O.UNIPROT),
            ("refseq", O.REFSEQ),
            ("rcsb pdb", O.PDBR),
            ("wwpdb", O.PDBR),
            ("taxid", O.TAXON),
        )
        for prefix, uri in PREFIXES:
            if not string.startswith(prefix + ":"):
                continue
            # remove the PSI-MI parens annotation and get the body
            string = string.split("(")[0]
            suffix = string[len(prefix) + 1:]
            # special cases
            if prefix == "psi-mi":
                suffix = suffix.replace("\"","").replace(":","_")
                return U(uri + suffix)
            elif prefix == "taxid":
                if "-" in suffix:
                    return L(suffix)
            elif prefix != "uniprotkb":
                uri += "_"
            # convert to URI
            return U(uri + suffix)
        raise ValueError("unknown PSI-MI prefix in '{}'".format(string))

    def _siphon_uris(self, triples, s, p, string, unique = False):
        """Converts a PSI-MI string list into triples.

        :param triples: the triples to add the results to.
        :param s: the triple subject.
        :param p: the triple predicate.
        :param string: the PSI-MI string list object(s).
        :param unique: checks that the string list has length 1.
        """
        parts = self._get_parts(string)
        if unique:
            assert len(parts) == 1
        triples.extend(map(lambda o: (s, p, self._prefix_to_uri(o)),
                           parts))

    def _siphon_confidences(self, triples, s, string):
        """Converts a PSI-MI `score` into triples.

        It only cares about the intact-miscore. Other scores are
        discarded.

        :param s: the interaction URI.
        :param string: the field value string.
        """
        triples.extend(map(lambda part: (s, O.OBO.MI_has_confidence, L(part[15:])),
                           filter(lambda part: part.startswith("intact-miscore:"),
                                  self._get_parts(string))))

    def _siphon_bool(self, triples, s, p, string,
                     skip_true = False, skip_false = False):
        o = { "true" : True, "false" : False }[string]
        if o == True and skip_true:
            return
        elif o == False and skip_false:
            return
        else:
            triples.append((s, p, L(o)))

    def _siphon_row(self, triples, row):
        """Converts a single PSI-MI row into RDF triples.

        Each interaction is reified by a blank node."""
        interaction = B()
        triples.append((interaction, O.RDF.type, O.OBO.MI_0000))

        OPS = {
            # PSI-MI TAB 2.5
            "ID_INTERACTOR_A" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.MI_has_interactor_a,
                                                value, unique = True),
            "ID_INTERACTOR_B" : 
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.MI_has_interactor_b,
                                                value, unique = True),
            "INTERACTION_DETECTION_METHOD" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.MI_has_detection_method,
                                                value),
            "TAXID_INTERACTOR_A" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.has_taxid_interactor_a, value),
            "TAXID_INTERACTOR_B" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.has_taxid_interactor_b, value),
            "INTERACTION_TYPES" :
                lambda value: self._siphon_uris(triples,
                                                interaction,
                                                O.OBO.MI_has_type,
                                                value),
            "SOURCE_DATABASE" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.MI_has_source_db,
                                                value),
            "INTERACTION_IDENTIFIERS" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.MI_has_identifier,
                                                value),
            "CONFIDENCE_VALUES" :
                lambda value: self._siphon_confidences(triples, interaction,
                                                       value),
            # PSI-MI TAB 2.6
            "TAXID_HOST" :
                lambda value: self._siphon_uris(triples, interaction,
                                                O.OBO.has_taxid_host, value),
            "IS_NEGATIVE" :
                lambda value: self._siphon_bool(triples, interaction,
                                                O.OBO.is_negative, value),
            # PSI-MI TAB 2.7
        }
        for field, func in OPS.items():
            func(row[field])

