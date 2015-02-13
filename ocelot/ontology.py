# -*- coding: utf-8 -*-

from rdflib import Namespace as NS, URIRef as U

#
# Namespaces
#

OCELOT  = NS("http://ocelot.disi.unitn.it/term#")
UNKNOWN = NS("http://ocelot.disi.unitn.it/unknown-term#")
RDF     = NS("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
RDFS    = NS("http://www.w3.org/2000/01/rdf-schema#")
OWL     = NS("http://www.w3.org/2002/07/owl#")
DC      = NS("http://purl.org/dc/elements/1.1/")
DCTERMS = NS("http://purl.org/dc/terms/")
GO      = NS("http://www.geneontology.org/go#")
OBO     = NS("http://purl.obolibrary.org/obo/")
MI      = NS("http://purl.obolibrary.org/obo/mi#")
REFSEQ  = NS("http://purl.obolibrary.org/obo/refseq#")
UNICORE = NS("http://purl.uniprot.org/core/")
UNIPROT = NS("http://purl.uniprot.org/uniprot/")
TAXON   = NS("http://purl.uniprot.org/taxonomy/")
PDBO    = NS("http://rdf.wwpdb.org/schema/pdbx-v40.owl#")
PDBR    = NS("http://rdf.wwpdb.org/pdb/")

BINDINGS = (
    ("ocelot",  OCELOT),
    ("rdf",     RDF),
    ("rdfs",    RDFS),
    ("owl",     OWL),
    ("dc",      DC),
    ("dcterms", DCTERMS),
    ("go",      GO),
    ("obo",     OBO),
    ("mi",      MI),
    ("unicore", UNICORE),
    ("uniprot", UNIPROT),
    ("taxon",   TAXON),
    ("pdbo",    PDBO),
    ("pdbr",    PDBR),
)

#
# URIs
#

# A GO ID, e.g. 
GO_ID                   = U(GO)

# A UniProt/SwissProt accession
UNIPROT_ID              = U(UNIPROT)

# A UniProt taxon
UNIPROT_TAXON           = U(TAXON)

# An SGD ID, e.g. `S000006060`
SGD_ID                  = OCELOT.sgd_id
SGD_ID_HAS_TYPE         = OCELOT.sgd_id_has_type
SGD_ID_QUALIFIER        = OCELOT.sgd_id_qualifier
SGD_ID_HAS_QUALIFIER    = OCELOT.sgd_id_has_qualifier
SGD_ID_IN_CHROMOSOME    = OCELOT.sgd_id_in_chromosome
SGD_ID_IN_STRAND        = OCELOT.sgd_id_in_strand
SGD_ID_STARTS_AT        = OCELOT.sgd_id_starts_at
SGD_ID_STOPS_AT         = OCELOT.sgd_id_stops_at
SGD_ID_HAS_XREF         = OCELOT.sgd_id_has_xref
SGD_ID_HAS_SEQUENCE     = OCELOT.sgd_id_has_sequence
SGD_ID_HAS_GOSLIM       = OCELOT.sgd_id_has_goslim

# An SGD feature name, e.g. `YAL001C`
SGD_FEATURE             = OCELOT.sgd_feature
SGD_FEATURE_TYPE        = OCELOT.sgd_feature_type

# An SGD chromosome, e.g. `12`
SGD_CHROMOSOME          = OCELOT.sgd_chromosome

# An individual SGD-PDB homology match
SGD_PDB_HOMOLOGY        = OCELOT.sgd_pdb_homology
SGD_PDB_ALIGNMENT       = OCELOT.sgd_pdb_alignment
SGD_PDB_SCORE           = OCELOT.sgd_pdb_score
SGD_PDB_HAS_QUERY       = OCELOT.sgd_pdb_has_query
SGD_PDB_QUERY_START     = OCELOT.sgd_pdb_query_start
SGD_PDB_QUERY_STOP      = OCELOT.sgd_pdb_query_stop
SGD_PDB_HAS_TARGET      = OCELOT.sgd_pdb_has_target
SGD_PDB_TARGET_START    = OCELOT.sgd_pdb_target_start
SGD_PDB_TARGET_STOP     = OCELOT.sgd_pdb_target_stop

# An individual SGD-InterPro domain hit
SGD_IPR_HIT             = OCELOT.sgd_ipr_hit
SGD_IPR_HIT_STARTS_AT   = OCELOT.sgd_ipr_hit_starts_at
SGD_IPR_HIT_STOPS_AT    = OCELOT.sgd_ipr_hit_stops_at
SGD_IPR_HIT_HAS_EVALUE  = OCELOT.sgd_ipr_hit_has_evalue
SGD_IPR_HIT_IS_TRUE     = OCELOT.sgd_ipr_hit_is_true

# An individual SGD interaction record
SGD_INT                 = OCELOT.sgd_int
SGD_INT_HAS_BAIT        = OCELOT.sgd_int_has_bait
SGD_INT_HAS_HIT         = OCELOT.sgd_int_has_hit
SGD_INT_TYPE            = OCELOT.sgd_int_type
SGD_INT_EXP_TYPE        = OCELOT.sgd_int_exp_type
SGD_INT_CUR_TYPE        = OCELOT.sgd_int_cur_type
SGD_INT_SOURCE          = OCELOT.sgd_int_source
SGD_INT_HAS_TYPE        = OCELOT.sgd_int_has_type
SGD_INT_HAS_EXP_TYPE    = OCELOT.sgd_int_has_exp_type
SGD_INT_HAS_CUR_TYPE    = OCELOT.sgd_int_has_cur_type
SGD_INT_HAS_SOURCE      = OCELOT.sgd_int_has_source

# The type of a PDB ID
PDBID                   = PDBO.PDBID
PDBID_CHAIN             = PDBO.PDBID_CHAIN
PDBID_CHAIN_HAS_PS      = OCELOT.pdb_id_has_ps
PDBID_CHAIN_HAS_SS      = OCELOT.pdb_id_has_ss
PDB_HAS_RESOLUTION      = OCELOT.pdb_has_resolution
PDB_HAS_EXP_METHOD      = OCELOT.pdb_has_exp_method

# DSSP secondary-structure type
DSSP_ALPHA_HELIX        = OCELOT.dssp_alpha_helix
DSSP_BETA_BRIDGE        = OCELOT.dssp_beta_bridge
DSSP_BETA_LADDER        = OCELOT.dssp_beta_ladder
DSSP_3_HELIX            = OCELOT.dssp_3_helix
DSSP_5_HELIX            = OCELOT.dssp_5_helix
DSSP_TURN               = OCELOT.dssp_turn
DSSP_BEND               = OCELOT.dssp_bend
DSSP_COIL               = OCELOT.dssp_coil

# A Pfam protein family
PFAM_ID                 = OCELOT.pfam_id
PFAM_TYPE               = OCELOT.pfam_type
PFAM_ID_HAS_TYPE        = OCELOT.pfam_id_has_type
PFAM_RESIDUE            = OCELOT.pfam_residue
PFAM_RESIDUE_OCCURS_IN  = OCELOT.pfam_residue_occurs_in
PFAM_RESIDUE_SEQNUM     = OCELOT.pfam_residue_seqnum
PFAM_RESIDUE_PS         = OCELOT.pfam_residue_ps
PFAM_RESIDUE_SS         = OCELOT.pfam_residue_ss

# An iPfam region
IPFAM_REGION            = OCELOT.ipfam_region
IPFAM_REGION_INSTANCE_OF= OCELOT.ipfam_region_instance_of
IPFAM_REGION_OCCURS_IN  = OCELOT.ipfam_region_occurs_in
IPFAM_REGION_STARTS_AT  = OCELOT.ipfam_region_starts_at
IPFAM_REGION_STOPS_AT   = OCELOT.ipfam_region_stops_at

# An iPfam family interaction
IPFAM_INT               = OCELOT.ipfam_int
IPFAM_INT_TYPE          = OCELOT.ipfam_int_type
IPFAM_INT_HAS_PFAM      = OCELOT.ipfam_int_has_pfam
IPFAM_INT_HAS_INT_TYPE  = OCELOT.ipfam_int_has_int_type

# An iPfam region interaction
IPFAM_REGION_INT            = OCELOT.ipfam_region_int
IPFAM_REGION_INT_OCCURS_IN  = OCELOT.ipfam_region_int_occurs_in
IPFAM_REGION_INT_HAS_REGION = OCELOT.ipfam_region_int_has_region

# SIFTS
SIFTS_SAME_AS           = OCELOT.sifts_same_as

# Generic 'occurs in' predicate
YIP_PROTEIN             = OCELOT.yip_protein
YIP_DOMAIN              = OCELOT.yip_domain
YIP_RESIDUE             = OCELOT.yip_residue
YIP_RESIDUE_HAS_POSITION= OCELOT.yip_residue_has_position
YIP_INTERACTS_WITH      = OCELOT.yip_interacts_with
YIP_NOT_INTERACTS_WITH  = OCELOT.yip_not_interacts_with
YIP_PARENT_OF           = OCELOT.yip_parent_of
YIP_INSTANCE_OF         = OCELOT.yip_instance_of

def uri(root, name):
    return U(root + u"." + name)

def pdb_id_to_uri(pdb_id, chain = None):
    if chain:
        pdb_id += "_" + chain
    return U(PDBR + pdb_id)

def go_to_uri(go_id):
    return U(GO_ID + go_id)

_DSSP_TO_URI = {
    "H" : DSSP_ALPHA_HELIX,
    "B" : DSSP_BETA_BRIDGE,
    "E" : DSSP_BETA_LADDER,
    "G" : DSSP_3_HELIX,
    "I" : DSSP_5_HELIX,
    "T" : DSSP_TURN,
    "S" : DSSP_BEND,
    "C" : DSSP_COIL,
    " " : DSSP_COIL,
}

def dssp_to_uri(code):
    return _DSSP_TO_URI[code]
