# -*- coding: utf8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter, iterate_csv
from ocelot.services import FASTA

from rdflib import URIRef as U, BNode as B, Literal as L
from glob import glob

class PDBConverter(Converter):
    """Converter for PDB (Protein Data Bank) [1] to RDF.
    
    Note that this converter does not actually convert the *structure* data,
    but only sequences and information about the database.

    *References*
    
    [1] http://www.rcsb.org/pdb/
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("ids",     self._siphon_ids),
            ("info",    self._siphon_info),
            ("ps-ss",   self._siphon_ps_ss),
        )
        super(PDBConverter, self).__init__(SUBTARGETS, *args, **kwargs)
        self.method = kwargs.get("method", "biopython")

    def _get_path(self, filename):
        import os
        return os.path.join(self.src, "PDB", filename)

    def _siphon_ids(self, triples):
        for path in glob(self._get_path("divided/mmCIF")):
            pdb_id = O.uri(O.PDB_ID, f.split(".")[0])
            triples.append((pdb_id, O.RDF.type, O.PDBID))

    def _siphon_info_mmlib(self, triples, path):
        pass

    def _siphon_info_biopython(self, triples, path):
        pass

#    def get_info_one_mmlib(self, path, triples):
#        from mmLib import FileIO as io
#
#        fp = io.OpenFile(path, "r")
#
#        pdb_id = path.split("/")[-1].split(".")[0]
#        for chain in io.LoadStructure(file = fp):
#            triples.extend([
#                (self.pdb_to_uri(pdb_id, chain.chain_id), NS.RDF.type, OSBR.PDBID_CHAIN.value),
#            ])
#
#    def get_info_one_biopython(self, path, triples):
#        from Bio import PDB
#        import gzip
#
#        with gzip.open(path, "rb") as ifp:
#            with open("temp.cif", "wb") as ofp:
#                ofp.write(ifp.read())
#
#        parser = PDB.MMCIFParser()
#        structure = parser.get_structure("whatever", "temp.cif")
#
#        pdb_id = path.split("/")[-1].split(".")[0]
#        for model in structure:
#            for chain in model:
#                triples.extend([
#                    (self.pdb_to_uri(pdb_id, chain.id), NS.RDF.type, OSBR.PDBID_CHAIN.value),
#                ])
#
#            # We only care about the first model; they supposedly
#            # all have the same number of chains...
#            break

    def _siphon_info(self, triples):
        if self.method == "biopython":
            siphon_one = self._siphon_info_biopython
        elif self.method == "mmlib":
            siphon_one = self._siphon_info_mmlib
        else:
            raise ValueError("invalid method '{}'".format(self.method))
        for path in glob(self._get_path("divided/mmCIF")):
            siphon_one(triples, path)

    def _siphon_ps_ss(self, triples):
        ps, ss = {}, {}
        for header, sequence in FASTA().read(self._get_path("ss.txt")):
            pdbid, chain, content = header.split(":")

            key = (pdbid, chain)
            if content == "sequence":
                assert not key in ps
                ps[key] = sequence
            elif content == "secstr":
                assert not key in ss
                ss[key] = sequence
            else:
                raise RuntimeError("invalid entry '{}, {}'".format(header, sequence))
        for (pdb_id, chain), sequence in ps.items():
            triples.append((self._pdb_uri(pdb_id, chain),
                            O.PDBID_CHAIN_HAS_PS, L(sequence)))
        for (pdb_id, chain), sequence in ss.items():
            triples.append((self._pdb_uri(pdb_id, chain),
                            O.PDBID_CHAIN_HAS_SS, L(sequence)))
