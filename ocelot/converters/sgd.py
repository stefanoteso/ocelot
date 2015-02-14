# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv
from ocelot.services import FASTA

from rdflib import URIRef as U, BNode as B, Literal as L

class SGDConverter(Converter):
    """Converter for the Saccharomyces Genome Database[1].
    
    The converter assumes that the `data` directory includes an `SGD` directory
    with the following tab-separated files taken directly from the SGD website:
    
      chromosome_length.tab         : not used
      dbxref.tab                    : DONE
      domains.tab                   : DONE
      gene_association.sgd          : not used
      go_protein_complex_slim.tab   : not used
      go_slim_mapping.tab           : DONE
      interaction_data.tab          : DONE
      pdb_homologs.tab              : DONE
      phenotype_data.tab            : not used
      protein_properties.tab        : not used
      SGD_CDS_xref.txt              : not used

    Please note that the interactions annotated in SGD are taken from BioGRID.
    
    *References*

    [1] http://www.yeastgenome.org/
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("features",        self._siphon_features),
            ("xref",            self._siphon_xrefs),
            ("sequences",       self._siphon_sequences),
            ("go-slim",         self._siphon_goslims),
            ("domains",         self._siphon_domains),
            ("pdb-homologues",  self._siphon_pdb_homologues),
            ("interactions",    self._siphon_interactions),
        )
        super(SGDConverter, self).__init__(SUBTARGETS, *args, **kwargs)

    def _get_path(self, basename):
        import os
        return os.path.join(self.src, "SGD", basename)

    @staticmethod
    def _sanitize(string):
        REPS = (
            (" ", "_"),
            ("'", "_SQUOTE_"),
            ('"', "_DQUOTE_"),
            ("(", "_"),
            (")", "_")
        )
        for f, t in REPS:
            string = string.replace(f, t)
        return string

    def _sanitize_all(self, string):
        return set(map(self._sanitize, string.split("|")))

    @staticmethod
    def _pdb_to_uri(string):
        if len(string) == 4:
            u, t = string.lower(), O.PDBO.pdb_id
        elif len(string) == 6:
            parts = string.split("_")
            u, t = parts[0].lower() + "_" + parts[1], O.PDBO.pdb_id_chain
        else:
            raise ValueError("unexpected string '{}'".format(string))
        return U(O.PDBR + u), t

    def _siphon_features(self, triples):
        """Converts the `SGD_features.tab` file."""
        # XXX there may be chromosomes not marked as features
        FIELDS = (
            "SGD_ID",           # Primary SGD ID
            "FEAT_TYPE",        # Feature type
            "FEAT_QUALIFIER",   # Feature qualifier (optional)
            "FEAT_NAME",        # Feature name (optional)
            "STD_GENE_NAME",    # Standard gene name (optional)
            "ALIASES",          # Aliases (optional, |-separated)
            "PARENT_FEAT_NAME", # Parent feature name (optional)
            "OTHER_SGD_IDS",    # Secondary SGD IDs (optional, |-separated)
            "CHROMOSOME",       # Chromosome (optional)
            "START_COORD",      # Start coordinate (optional)
            "STOP_COORD",       # Stop coordinate (optional)
            "STRAND",           # Strand (optional)
            "GENETIC_POS",      # Genetic position (optional)
            "COORD_VERSION",    # Coordinate version (optional)
            "SEQ_VERSION",      # Sequence version (optional)
            "DESCRIPTION",      # Description (optional)
        )

        feat_types, qualifiers, chromosomes = set(), set(), set()
        for row in iterate_csv(self._get_path("SGD_features.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            sgd_id      = O.uri(O.SGD_ID, self._sanitize(row["SGD_ID"]))
            feat_id     = O.uri(O.SGD_FEATURE, self._sanitize(row["FEAT_NAME"]))
            feat_type   = O.uri(O.SGD_FEATURE_TYPE, self._sanitize(row["FEAT_TYPE"]))
            triples.extend([
                (sgd_id,    O.RDF.type,         O.SGD_ID),
                (feat_id,   O.RDF.type,         O.SGD_FEATURE),
                (sgd_id,    O.OWL.sameAs,       feat_id),
                (sgd_id,    O.SGD_ID_HAS_TYPE,  feat_type),
            ])
            self._add_if_new(triples, feat_types, feat_type,
                             O.RDF.type, O.SGD_FEATURE_TYPE)
            for q in filter(lambda q: len(q), self._sanitize_all(row["FEAT_QUALIFIER"])):
                q = O.uri(O.SGD_ID_QUALIFIER, q)
                self._add_if_new(triples, qualifiers, q,
                                 O.RDF.type, O.SGD_ID_QUALIFIER)
                triples.append((sgd_id, O.SGD_ID_HAS_QUALIFIER, q))
            for a in filter(lambda a: len(a), self._sanitize_all(row["OTHER_SGD_IDS"])):
                triples.append((sgd_id, O.OWL.sameAs, O.uri(O.SGD_ID, a)))
            parent = row["PARENT_FEAT_NAME"]
            if len(parent):
                parent = O.uri(O.SGD_FEATURE, self._sanitize(parent))
                triples.extend([
                    (parent,    O.RDF.type,         O.SGD_FEATURE),
                    (sgd_id,    O.DCTERMS.isPartOf, parent),
                ])
            chromosome = row["CHROMOSOME"]
            if len(chromosome):
                chromosome = O.uri(O.SGD_FEATURE, "chrosomosome_" + chromosome)
                self._add_if_new(triples, chromosomes, chromosome,
                                 O.RDF.type, O.SGD_CHROMOSOME)
                triples.append((sgd_id, O.SGD_ID_IN_CHROMOSOME, chromosome))
            if len(row["STRAND"]):
                triples.append((sgd_id, O.SGD_ID_IN_STRAND, L(row["STRAND"])))
            if len(row["START_COORD"]):
                triples.append((sgd_id, O.SGD_ID_STARTS_AT, L(int(row["START_COORD"]))))
            if len(row["STOP_COORD"]):
                triples.append((sgd_id, O.SGD_ID_STOPS_AT, L(int(row["STOP_COORD"]))))

    def _siphon_xrefs(self, triples):
        """Converts the `dbxref.tab` file."""
        FIELDS = (
            "XREF_ID",          # Cross-reference ID
            "XREF_ID_SOURCE",   # Cross-reference database ID
            "XREF_ID_TYPE",     # Cross-reference type (like 'PDB chain' or 'PDB best hit')
            "FEAT_NAME",        # ORF name
            "SGD_ID",           # SGD ID
            "UNDOCUMENTED",     # Undocumented
        )
        for row in iterate_csv(self._get_path("dbxref.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            sgd_id      = O.uri(O.SGD_ID, self._sanitize(row["SGD_ID"]))
            xref_source = self._sanitize(row["XREF_ID_SOURCE"])
            xref_type   = self._sanitize(row["XREF_ID_TYPE"])
            if (xref_source, xref_type) == ("EBI", self._sanitize("UniProt/TrEMBL ID")) or \
               (xref_source, xref_type) == ("EBI", self._sanitize("UniProt/Swiss-Prot ID")):
                xref_id = O.UNIPROT_ID + row["XREF_ID"]
            else:
                continue
            triples.append((sgd_id, O.SGD_ID_HAS_XREF, xref_id))

    def _siphon_sequences(self, triples):
        """Converts the `orf_trans_all.fasta` file."""
        for header, sequence in FASTA().read(self._get_path("orf_trans_all.fasta")):
            assert "SGDID:" in header

            sgd_ids = filter(lambda word: word.startswith("SGDID:"),
                             header.split())
            assert len(sgd_ids) == 1

            sgd_id = O.uri(O.SGD_ID, self._sanitize(sgd_ids[0].split(":")[1].strip(",")))
            triples.append((sgd_id, O.SGD_ID_HAS_SEQUENCE, L(sequence)))

    def _siphon_domains(self, triples):
        """Converts the `domains.tab` file.

        The data comes from an InterPro scan over the SGD entries.
        """
        FIELDS = (
            "FEAT_NAME",        # S. cerevisiae systematic name (ID of the input sequence)
            "CRC64",            # CRC of the proteic sequence
            "LENGTH",           # Lenght of the sequence in AA
            "METHOD",           # Analysis method
            "DB_MEMBERS",       # DB members entry for this match
            "DB_DESCRIPTION",   # DB member description for the entry
            "START",            # start of the domain match
            "STOP",             # end of the domain match
            "EVALUE",           # E-value of the match (defined by DB)
            "STATUS",           # Status of the match: T=true, ?=unknown
            "DATE",             # Date of the run
            "IPR_ID",           # InterPro ID
            "IPR_DESCRIPTION",  # InterPro description
            "IPR_GO",           # GO description of the InterPro entry
        )
        for row in iterate_csv(self._get_path("domains.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            feat_id = O.uri(O.SGD_FEATURE, self._sanitize(row["FEAT_NAME"]))
            is_true = L({ "T":True, "?":False }[row["STATUS"]])

            try:
                evalue = L(float(row["EVALUE"]))
            except ValueError:
                evalue = L(-1.0)

            _ = B()
            triples.extend([
                (_, O.RDF.type,                 O.SGD_IPR_HIT),
                (_, O.SGD_IPR_HIT_STARTS_AT,    L(int(row["START"]))),
                (_, O.SGD_IPR_HIT_STOPS_AT,     L(int(row["STOP"]))),
                (_, O.SGD_IPR_HIT_HAS_EVALUE,   evalue),
                (_, O.SGD_IPR_HIT_IS_TRUE,      is_true),
            ])

    def _siphon_goslims(self, triples):
        """Converts the `go_slim_mapping.tab` file."""
        # XXX there are a few goslim annotations with commas
        FIELDS = (
            "GENE_ORF"      # Systematic gene name
            "GENE_STD",     # Gene name (optional)
            "SGD_ID",       # Gene SGD ID
            "GO_ASPECT",    # P = Process, F = Function, C = Component
            "GO_TERM",      # GO SLIM term
            "GO_ID",        # GO term ID
            "FEATURE_TYPE", # Such as 'ORF' or 'tRNA'
        )
        for row in iterate_csv(self._get_path("go_slim_mapping.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            sgd_id  = O.uri(O.SGD_ID, self._sanitize(row["SGD_ID"]))
            goslim  = O.go_to_uri(self._sanitize(row["GO_ID"]))
            triples.append((sgd_id, O.SGD_ID_HAS_GOSLIM, goslim))

    def _siphon_pdb_homologues(self, triples):
        """Converts the `pdb_homologs.tab` file."""
        FIELDS = (
            "FEAT_NAME",            # S. cerevisiae systematic name
            "START_COORD_QUERY",    # start coord (aa position) in yeast
            "STOP_COORD_QUERY",     # stop coord (aa position) in yeast
            "START_COORD_TARGET",   # start coord (aa position) in target
            "STOP_COORD_TARGET",    # stop coord (aa position) in target
            "PERCENT_ALIGNED",      # percent of yeast contained in target
            "SCORE",                # log of the expectation value
            "TARGET_PDB_ID",        # PDB identifier
            "TARGET_TAXON_ID",      # target taxon ID
            "TARGET_TAXON_NAME",    # target taxon species name
        )
        for row in iterate_csv(self._get_path("pdb_homologs.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            pdb_id_chain, _ = self._pdb_to_uri(row["TARGET_PDB_ID"])

            match = B()
            triples.extend([
                (match,     O.RDF.type,             O.SGD_PDB_HOMOLOGY),
                (match,     O.SGD_PDB_HAS_QUERY,    O.uri(O.SGD_FEATURE, row["FEAT_NAME"])),
                (match,     O.SGD_PDB_ALIGNMENT,    L(float(row["PERCENT_ALIGNED"]) / 100.0)),
                (match,     O.SGD_PDB_HAS_TARGET,   pdb_id_chain),
            ])

    def _siphon_interactions(self, triples):
        """Converts the `interaction_data.tab` file."""
        FIELDS = (
            "BAIT_FEAT_NAME",   # Bait: feature (ORF) name and gene name (optional)
            "BAIT_STD_NAME",    #
            "HIT_FEAT_NAME",    # Hit: feature (ORF) name and gene name (optional)
            "HIT_STD_NAME",     #
            "EXPERIMENT_TYPE",  # Description of the experiment
            "INTERACTION_TYPE", # 'Genetic' or 'Physical'
            "SOURCE_DATABASE",  # Source database
            "CURATION_TYPE",    # Manual or high-throughput
            "NOTES",            # Free text (useless)
            "PHENOTYPE",        # Phenotype of the interaction (optional)
            "REFERENCE",        # List of references as 'SGD_REF:" (SGDID) or 'PMID:' (PubMed)
            "CITATION",         # List of citations
        )
        int_types, exp_types, cur_types, sources = set(), set(), set(), set()
        for row in iterate_csv(self._get_path("interaction_data.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            bait        = O.uri(O.SGD_FEATURE, self._sanitize(row["BAIT_FEAT_NAME"]))
            hit         = O.uri(O.SGD_FEATURE, self._sanitize(row["HIT_FEAT_NAME"]))
            int_type    = O.uri(O.SGD_INT_TYPE, self._sanitize(row["INTERACTION_TYPE"]))
            exp_type    = O.uri(O.SGD_INT_EXP_TYPE, self._sanitize(row["EXPERIMENT_TYPE"]))
            cur_type    = O.uri(O.SGD_INT_CUR_TYPE, self._sanitize(row["CURATION_TYPE"]))
            source      = O.uri(O.SGD_INT_SOURCE, self._sanitize(row["SOURCE_DATABASE"]))

            self._add_if_new(triples, int_types, int_type, O.RDF.type, O.SGD_INT_TYPE)
            self._add_if_new(triples, exp_types, exp_type, O.RDF.type, O.SGD_INT_EXP_TYPE)
            self._add_if_new(triples, cur_types, cur_type, O.RDF.type, O.SGD_INT_CUR_TYPE)
            self._add_if_new(triples, source, source, O.RDF.type, O.SGD_INT_SOURCE)

            interaction = B()
            triples.extend([
                (interaction, O.RDF.type, O.SGD_INT),
                (interaction, O.SGD_INT_HAS_BAIT, bait),
                (interaction, O.SGD_INT_HAS_HIT, hit),
                (interaction, O.SGD_INT_HAS_TYPE, int_type),
                (interaction, O.SGD_INT_HAS_EXP_TYPE, exp_type),
                (interaction, O.SGD_INT_HAS_CUR_TYPE, cur_type),
                (interaction, O.SGD_INT_HAS_SOURCE, source),
            ])
