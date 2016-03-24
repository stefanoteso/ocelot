# -*- coding: utf-8 -*-

import os

import ocelot.ontology as O
from ocelot.converters.base import Converter
from ocelot.sequences import read_fasta
from ocelot.services import iterate_csv

from rdflib import URIRef as U, Literal as L, BNode as B

class CAFA13Converter(Converter):
    """Converter for the CAFA13 dataset.

    .. todo::

        Add support for the *CyC annotations.
        Add support for the test data.

    The converter assumes that the ``CAFA13`` directory includes WRITEME.
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("go",        self._siphon_go),
            ("swissprot", self._siphon_swissprot),
            ("uniprot",   self._siphon_uniprot),
        )
        super(CAFA13Converter, self).__init__(SUBTARGETS, *args, **kwargs)

    def _join(self, *args):
        return os.path.join(self.src, "CAFA13", "CAFATrainingData", *args)

    def _siphon_go(self, triples):
        # p_seq = read_fasta(path))
        # assert len(p_seq) = 118411
        # for p, seq in p_seq:
        #     p = O.uri(O.CAFA13_ID, p)
        #     triples.extend([
        #         (p, O.RDF.type, O.CAFA13_ID),
        #         (p, O.CAFA13_ID_HAS_SEQ, L(seq)),
        #     ])

        # FIELDS = ("UNIPROT_ID", "TAXON", "ANNOTATION")
        # num_rows = 0
        # for row in iterate_csv(self._join("GO", "go_function.dat"),
        #                        delimiter = "\t", fieldnames = FIELDS):
        #     uniprot_id  = U(O.UNIPROT_ID + row["UNIPROT_ID"])
        #     taxon       = U(O.TAXON + row["TAXON"])
        #     term_id, namespace, evidence_code = row["ANNOTATION"].split(",")
        #     annotation = B()
        #     triples.extend([
        #         (uniprot_id, O.RDF.type, O.CAFA13_ID),
        #         (uniprot_id, O.CAFA13_ID_HAS_ANNOT, annotation),
        #         (annotation, O.RDF.type, O.CAFA13_ANNOT),
        #         (annotation, O.CAFA13_ANNOT_HAS_GO, U(O.GO_ID + term_id)),
        #         (annotation, O.CAFA13_ANNOT_HAS_TAXID, taxon),
        #         (annotation, O.CAFA13_ANNOT_HAS_EVCODE, L(evidence_code)),
        #     ])
        #     num_rows += 1
        # assert num_rows == 760933
        pass

    def _siphon_swissprot(self, triples):
        # FIELDS = ("UNIPROT_NAME", "TAXON", "UNIPROT_IDS", "ANNOTATIONS")
        # num_rows = 0
        # for row in iterate_csv(self._join("GO", "go_function.dat"),
        #                        delimiter = "\t", fieldnames = FIELDS):
        #     name        = L(row["UNIPROT_NAME")
        #     taxon       = U(O.TAXON + row["TAXON"])
        #     annotation = B()
        #     for part in row["UNIPROT_IDS"].split(";"):
        #         uniprot_id = U(O.UNIPROT_ID + part)
        #         triples.extend([
        #             (uniprot_id, O.RDF.type, O.CAFA13_ID),
        #             (uniprot_id, O.CAFA13_ID_HAS_ANNOT, annotation),
        #         ])
        #     for part in row["ANNOTATIONS"].split(";"):
        #         term_id, namespace, evidence_code = part.split(",")
        #             (annotation, O.RDF.type, O.CAFA13_ANNOT),
        #             (annotation, O.CAFA13_ANNOT_HAS_GO, U(O.GO_ID + term_id)),
        #             (annotation, O.CAFA13_ANNOT_HAS_TAXID, taxon),
        #             (annotation, O.CAFA13_ANNOT_HAS_EVCODE, L(evidence_code)),
        #         ])
        #     num_rows += 1
        # assert num_rows == 55061
        pass

    def _siphon_uniprot(self, triples):
        # # EVIDENCE CODES
        # # --------------
        # #
        # #  EXP: Inferred from Experiment
        # #  IC: Inferred by Curator
        # #  IDA: Inferred from Direct Assay
        # #  IEA: Inferred from Electronic Annotation
        # #  IEP: Inferred from Expression Pattern
        # #  IGC: Inferred from Genomic Context
        # #  IGI: Inferred from Genetic Interaction
        # #  IMP: Inferred from Mutant Phenotype
        # #  IPI: Inferred from Physical Interaction
        # #  ISS: Inferred from Sequence or Structural Similarity
        # #  NAS: Non-traceable Author Statement
        # #  RCA: inferred from Reviewed Computational Analysis
        # #  TAS: Traceable Author Statement
        # #
        # # See `this <http://www.uniprot.org/help/gene_ontology>`_.
        # FIELDS = ("UNIPROT_ID", "TAXON", "QUALIFIER", "ANNOTATION")
        # QUALIFIERS = ["", "NOT",
        #               "colocalizes_with", "NOT|colocalizes_with",
        #               "contributes_to", "NOT|contributes_to"]
        # num_rows = 0
        # for row in iterate_csv(self._join("UniProt_GOA", "uniprot_goa_function.dat"),
        #                        delimiter = "\t", fieldnames = FIELDS):
        #     uniprot_id  = U(O.UNIPROT_ID + row["UNIPROT_ID"])
        #     taxon       = U(O.TAXON + row["TAXON"])
        #     qualifier   = row["QUALIFIER"]
        #     assert qualifier in QUALIFIERS
        #     term_id, namespace, evidence_code = row["ANNOTATION"].split(",")
        #     annotation = B()
        #     triples.extend([
        #         (uniprot_id, O.RDF.type, O.CAFA13_ID),
        #         (uniprot_id, O.CAFA13_ID_HAS_ANNOT, annotation),
        #         (annotation, O.RDF.type, O.CAFA13_ANNOT),
        #         (annotation, O.CAFA13_ANNOT_HAS_GO, U(O.GO_ID + term_id)),
        #         (annotation, O.CAFA13_ANNOT_HAS_TAXID, taxon),
        #         (annotation, O.CAFA13_ANNOT_HAS_QUAL, L(qualifier)),
        #         (annotation, O.CAFA13_ANNOT_HAS_EVCODE, L(evidence_code)),
        #     ])
        #     num_rows += 1
        # assert num_rows == 810380
        pass
