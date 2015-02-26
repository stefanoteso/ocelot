# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv

class SIFTSConverter(Converter):
    """Converter for `SIFTS <http://www.ebi.ac.uk/pdbe/docs/sifts/>`_.
    
    SIFTS provides mappings between IDs of different biological databases,
    including PDB, Uniprot, GO, Pfam, InterPro, SCOP and CATH.
    
    Tested with SIFTS from 2014-07-08.
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETs = (
            ("pdb-uniprot", self._siphon_pdb_uniprot)
            #("pdb-pfam",   self._siphon_pdb_pfam)
            #("pdb-taxon",  self._siphon_pdb_taxon)
        )
        super(SIFTSConverter, self).__init__(SUBTARGETS,
                                             *args, **kwargs)

    def _get_path(self, basename):
        import os
        return os.path.join(self.src, "SGD", basename)

    def _siphon_pdb_uniprot(self, triples):
        """Converts the `pdb_chain_uniprot.csv` file."""

        FIELDS = (
            "PDB_ID",
            "PDB_CHAIN",
            "SP_ID",
            "RES_BEGIN",
            "RES_END",
            "PDB_BEGIN",
            "PDB_END",
            "SP_BEGIN",
            "SP_END"
        )

        for row in iterate_csv(self._get_path("pdb_chain_uniprot.csv"),
                               delimiter = "\t", fieldnames = FIELDS):
            pdb_id_chain = self._pdb_uri(row["PDB_ID"], row["PDB_CHAIN"])
            sp_id = self._uniprot_uri(row["SP_ID"])

            triples.append((pdb_id_chain, O.SIFTS_SAME_AS, sp_id))
