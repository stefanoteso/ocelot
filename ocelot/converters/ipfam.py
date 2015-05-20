# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv

from rdflib import URIRef as U, BNode as B, Literal as L

class IPfamConverter(Converter):
    """Converter for `iPfam <http://www.ipfam.org/>`_.

    The converter assumes that the data directory includes an ``ipfam``
    directory with a full iPfam dump, laid out like this::

      homodomain_interactions.csv   : DONE
      heterodomain_interactions.csv : DONE
      protein_family.txt            : DONE
      pdb_entry.txt                 : DONE
      pdb_residue_data.txt          : DONE
      pdb_protein_region.txt        : DONE
      pdb_protein_region_int.txt    : DONE
      pdb_protein_res_int.txt       : DONE
      pdb_protein_atom_int.txt      : not used
      pdb_protein_region_lig_int.txt: not used
      pdb_protein_res_lig_int.txt   : not used
      pdb_protein_lig_atom_int.txt  : not used
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("pdb-entries",             self._siphon_pdb_entries),
            ("families",                self._siphon_families),
            ("family-interactions",     self._siphon_family_interactions),
            ("regions",                 self._siphon_regions),
            ("region-interactions",     self._siphon_region_interactions),
            ("residues",                self._siphon_residues),
            ("residue-interactions",    self._siphon_residue_interactions),
        )
        super(IPfamConverter, self).__init__(SUBTARGETS,
                                             *args, **kwargs)

    def _get_path(self, filename, is_db = True):
        import os
        parts = [ self.src, "ipfam", "1.0" ]
        if is_db:
            parts += [ "database_files" ]
        parts += [ filename ]
        return os.path.join(*parts)

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

    def _siphon_pdb_entries(self, triples):
        """Converts the `pdb_entry.txt` file."""

        FIELDS = (
            "PDB_ID",       # `pdb_id` varchar(4) NOT NULL DEFAULT 'NULL',
            "HEADER",       # `header` text,
            "TITLE",        # `title` text,
            "DATE",         # `date` date NOT NULL,
            "RESOLUTION",   # `resolution` decimal(5,2) unsigned NOT NULL,
            "EXP_METHOD",   # `expt_method` text NOT NULL,
            "AUTHOR",       # `author` mediumtext,
            "PDB_FILE",     # `pdb_file` int(10) DEFAULT '0',
            "SIFTS_FILE"    # `sifts_file` int(10) DEFAULT '0',
        )

        for row in iterate_csv(self._get_path("pdb_entry.txt"),
                               delimiter = "\t", fieldnames = FIELDS):
            pdb_id  = U(O.PDBR + row["PDB_ID"].lower())
            res     = L(float(row["RESOLUTION"]))
            method  = L(self._sanitize(row["EXP_METHOD"]))

            triples.extend([
                (pdb_id, O.RDF.type, O.PDBID),
                (pdb_id, O.PDB_HAS_RESOLUTION, res),
                (pdb_id, O.PDB_HAS_EXP_METHOD, method),
            ])

    def _siphon_families(self, triples):
        """Converts the `protein_family.txt file."""
        FIELDS = (
            "FAMILY_INT",       # `auto_prot_fam` int(11) NOT NULL AUTO_INCREMENT,
            "FAMILY_ACC",       # `accession` varchar(45) DEFAULT NULL,
            "FAMILY_ID",        # `identifier` varchar(45) DEFAULT NULL,
            "DESCRIPTION",      # `description` text,
            "COMMENT",          # `comment` longtext,
            "FAMILY_TYPE",      # `type` enum('family','domain','motif','repeat') DEFAULT NULL,
            "SOURCE_DB",        # `source_db` enum('pfama') DEFAULT NULL,
            "COLOUR",           # `colour` varchar(7) DEFAULT NULL,
            "NUMBER_FAM_INT",   # `number_fam_int` int(5) DEFAULT '0',
            "NUMBER_LIG_INT",   # `number_lig_int` int(5) DEFAULT '0',
            "NUMBER_PDBS",      # `number_pdbs` int(5) DEFAULT '0',
        )

        pfam_types = set()
        for row in iterate_csv(self._get_path("protein_family.txt"),
                               num_skip = 1,
                               delimiter = "\t", fieldnames = FIELDS):
            assert len(row["FAMILY_ACC"]) > 0

            if row["SOURCE_DB"] != "pfama":
                print "Warning: invalid row '{}'".format(row.items())
                continue

            pfam_acc = O.uri(O.PFAM_ID, row["FAMILY_ACC"])
            triples.append((pfam_acc, O.RDF.type, O.PFAM_ID))
            if row["FAMILY_TYPE"] is None:
                continue
            if row["FAMILY_TYPE"]:
                pfam_type = O.uri(O.PFAM_TYPE, row["FAMILY_TYPE"])
                triples.append((pfam_acc, O.PFAM_ID_HAS_TYPE, pfam_type))
                if not pfam_type in pfam_types:
                    pfam_types.add(pfam_type)
                    triples.append((pfam_type, O.RDF.type, O.PFAM_TYPE))

    def _siphon_family_interactions(self, triples):
        ints = []
        FIELDS = ("FAM", "ID", "INT_TYPE")
        for row in iterate_csv(self._get_path("homodomain_interaction.csv", is_db = False),
                               num_skip = 1, delimiter = "\t", fieldnames = FIELDS):
            ints.append((row["FAM"], row["FAM"], row["INT_TYPE"]))
        FIELDS = ("FAM1", "ID1", "FAM2", "ID2", "INT_TYPE")
        for row in iterate_csv(self._get_path("heterodomain_interaction.csv", is_db = False),
                               num_skip = 1, delimiter = "\t", fieldnames = FIELDS):
            ints.append((row["FAM1"], row["FAM2"], row["INT_TYPE"]))
        int_types = set()
        for fam1, fam2, int_type in ints:
            fam1 = O.uri(O.PFAM_ID, fam1)
            fam2 = O.uri(O.PFAM_ID, fam2)
            int_type = O.uri(O.IPFAM_INT_TYPE, int_type)
            blank = B()
            triples.extend([
                (blank, O.RDF.type,                 O.IPFAM_INT),
                (blank, O.IPFAM_INT_HAS_INT_TYPE,   int_type),
                (blank, O.IPFAM_INT_HAS_PFAM,       fam1),
                (blank, O.IPFAM_INT_HAS_PFAM,       fam2),
            ])
            self._add_if_new(triples, int_types,
                             int_type, O.RDF.type, O.IPFAM_INT_TYPE)
        del ints

    def _siphon_regions(self, triples):
        FIELDS = (
            "REGION",           # int(10)
            "PROT_FAM",         # int(11)
            "PROT_FAM_ACC",     # varchar(45)
            "PDB_ID",           # varchar(4)
            "CHAIN",            # varchar(1)
            "START",            # int(11)
            "START_ICODE",      # varchar(1)
            "END",              # int(11)
            "END_ICODE",        # varchar(1)
            "REGION_SOURCE_DB", # varchar(12)
        )
        for row in iterate_csv(self._get_path("pdb_protein_region.txt"),
                               delimiter = "\t", fieldnames = FIELDS):
            region          = O.uri(O.IPFAM_REGION, row["REGION"])
            pfam            = O.uri(O.PFAM_ID, row["PROT_FAM_ACC"])
            pdb_id_chain    = U(O.PDBR + row["PDB_ID"].lower() + "_" + row["CHAIN"])
            triples.extend([
                (region,    O.RDF.type, O.IPFAM_REGION),
                (region,    O.IPFAM_REGION_INSTANCE_OF, pfam),
                (region,    O.IPFAM_REGION_OCCURS_IN,   pdb_id_chain),
                (region,    O.IPFAM_REGION_STARTS_AT,   L(int(row["START"]))),
                (region,    O.IPFAM_REGION_STOPS_AT,    L(int(row["END"]))),
            ])

    def _siphon_region_interactions(self, triples):
        FIELDS = (
            "REGION_INT",       # `auto_reg_int` bigint(20) NOT NULL AUTO_INCREMENT,
            "PDB_ID",           # `pdb_id` varchar(4) NOT NULL,
            "REGION_A",         # `region_id_A` int(10) unsigned NOT NULL,
            "REGION_B",         # `region_id_B` int(10) unsigned NOT NULL,
            "IS_INTRACHAIN",    # `intrachain` tinyint(1) NOT NULL,
            "QUALITY_CONTROL",  # `quality_control` int(10) unsigned NOT NULL,
        )
        for row in iterate_csv(self._get_path("pdb_protein_region_int.txt"),
                               delimiter = "\t", fieldnames = FIELDS):
            region_int  = O.uri(O.IPFAM_REGION_INT, row["REGION_INT"])
            region_a    = O.uri(O.IPFAM_REGION, row["REGION_A"])
            region_b    = O.uri(O.IPFAM_REGION, row["REGION_B"])
            pdb_id      = U(O.PDBR + row["PDB_ID"].lower())
            triples.extend([
                (region_int,    O.RDF.type,                     O.IPFAM_REGION_INT),
                (region_int,    O.IPFAM_REGION_INT_OCCURS_IN,   pdb_id),
                (region_int,    O.IPFAM_REGION_INT_HAS_REGION,  region_a),
                (region_int,    O.IPFAM_REGION_INT_HAS_REGION,  region_b),
            ])

    def _siphon_residues(self, triples):
        FIELDS = (
            "PDB_ID",               # `pdb_id` varchar(4) NOT NULL,
            "CHAIN",                # `chain` char(1) NOT NULL,
            "SERIAL",               # `serial` int(10) DEFAULT NULL,
            "PDB_RES",              # `pdb_res` char(3) DEFAULT NULL,
            "PDB_SEQ_NUMBER",       # `pdb_seq_number` int(10) DEFAULT NULL,
            "PDB_INS_CODE",         # `pdb_insert_code` varchar(1) DEFAULT NULL,
            "OBSERVED",             # `observed` int(1) DEFAULT NULL,
            "DSSP_CODE",            # `dssp_code` varchar(1) DEFAULT NULL,
            "UNIPROT_ACC",          # `uniprot_acc` varchar(6) DEFAULT NULL,
            "UNIPROT_RES",          # `uniprot_res` char(3) DEFAULT NULL,
            "UNIPROT_SEQ_NUMBER",   # `uniprot_seq_number` int(10) DEFAULT NULL,
        )
        pass

    def _siphon_residue_interactions(self, triples):
        FIELDS = (
            "RESIDUE_INT",  # auto_res_int` bigint(20) NOT NULL AUTO_INCREMENT,
            "REGION_INT",   # auto_reg_int` bigint(20) NOT NULL,
            "PDB_ID",       # pdb_id` varchar(4) NOT NULL,
            "RESIDUE_A",    # residue_A` int(10) NOT NULL,
            "RESIDUE_B",    # residue_B` int(10) NOT NULL,
            "BOND_TYPE",    # bond` enum('vdw','hbond','hbback','hbside','electro','covalent') NOT NULL,
        )
        pass
