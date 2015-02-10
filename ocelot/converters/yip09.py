# -*- coding: utf-8 -*-

import ocelot.ontology as O
from ocelot.converters.base import Converter

from rdflib import URIRef as U, BNode as B, Literal as L

import os

# TODO: simplify a lot, e.g. read the files only once...
# TODO: check that the number of objects is correct

class Yip09Converter(Converter):
    """Converter for the Yip et al.[1] dataset.

    The converter assumes that the `data` directory includes the `ppi`, `ddi`
    and `rri` directories from the dataset hosted at [2].

    *References*

    [1] Yip, Kim, McDermott, Gerstein, "Multi-level learning: improving the
    prediction of protein, domain and residue interactions by allowing
    information flow between levels", BMC Bioinformatics, 2009.
    [2] http://networks.gersteinlab.org/mll
    """
    def __init__(self, *args, **kwargs):
        SUBTARGETS = (
            ("pin",         self._siphon_pin),
            ("din",         self._siphon_din),
            ("rin",         self._siphon_rin),
            ("parentpd",    self._siphon_parentpd),
            ("parentdr",    self._siphon_parentdr),
            ("aliases",     self._siphon_aliases),
        )
        super(Yip09Converter, self).__init__("yip", SUBTARGETS, *args, **kwargs)

    def _read(self, path):
        with open(os.path.join(self.src, "yip09", path)) as fp:
            return [ line.strip().split() for line in fp ]

    @staticmethod
    def _p_uri(p):
        return O.uri(O.YIP_PROTEIN, p)

    @staticmethod
    def _sgd_uri(p):
        return O.uri(O.SGD_FEATURE, p)

    @staticmethod
    def _d_uri(p, d):
        return O.uri(O.YIP_DOMAIN, p + "_" + d)

    @staticmethod
    def _pfam_uri(d):
        return O.uri(O.PFAM_ID, d)

    @staticmethod
    def _r_uri(p, d, r):
        return O.uri(O.YIP_RESIDUE, p + "_" + d + "_" + r)

    def _siphon_pin(self, triples):
        for words in self._read(os.path.join("ppi", "ready", "proteins.txt")):
            assert len(words) == 1
            triples.extend([
                (self._p_uri(words[0]), O.RDF.type, O.YIP_PROTEIN),
                (self._p_uri(words[0]), O.OWL.sameAs, self._sgd_uri(words[0])),
            ])
        pos_triples = [
            (self._p_uri(words[0]), O.YIP_INTERACTS_WITH, self._p_uri(words[1]))
            for words in self._read(os.path.join("ppi", "ready", "goldPosProteinPairs.txt"))]
        neg_triples = [
            (self._p_uri(words[0]), O.YIP_NOT_INTERACTS_WITH, self._p_uri(words[1]))
            for words in self._read(os.path.join("ppi", "ready", "goldNegProteinPairs.txt"))]
        assert len(pos_triples) == 3201
        triples.extend(pos_triples + neg_triples)

    def _siphon_din(self, triples):
        for words in self._read(os.path.join("ddi", "ready", "domains.txt")):
            assert len(words) == 2
            triples.extend([
                (self._d_uri(*words[0:2]), O.RDF.type, O.YIP_DOMAIN),
                (self._d_uri(*words[0:2]), O.YIP_INSTANCE_OF, self._pfam_uri(words[1])),
            ])
        for words in self._read(os.path.join("ddi", "ready", "goldPosDomainPairs.txt")):
            assert len(words) == 5
            triples.extend([
                (self._d_uri(*words[0:2]), O.YIP_INTERACTS_WITH, self._d_uri(*words[2:4])),
            ])
        for words in self._read(os.path.join("ddi", "ready", "goldNegDomainPairs.txt")):
            assert len(words) == 5
            triples.extend([
                (self._d_uri(*words[0:2]), O.YIP_NOT_INTERACTS_WITH, self._d_uri(*words[2:4])),
            ])

    def _siphon_rin(self, triples):
        for words in self._read(os.path.join("rri", "ready", "residues.txt")):
            assert len(words) == 3
            triples.extend([
                (self._r_uri(*words[0:3]), O.RDF.type, O.YIP_RESIDUE),
            ])
        for words in self._read(os.path.join("rri", "ready", "goldPosResiduePairs.txt")):
            assert len(words) == 7
            triples.extend([
                (self._r_uri(*words[0:3]), O.YIP_INTERACTS_WITH, self._r_uri(*words[3:6])),
            ])
        for words in self._read(os.path.join("rri", "ready", "goldNegResiduePairs.txt")):
            assert len(words) == 7
            triples.extend([
                (self._r_uri(*words[0:3]), O.YIP_NOT_INTERACTS_WITH, self._r_uri(*words[3:6])),
            ])

    def _siphon_parentpd(self, triples):
        for words in self._read(os.path.join("ddi", "ready", "domains.txt")):
            assert len(words) == 2
            triples.extend([
                    (self._p_uri(words[0]), O.YIP_PARENT_OF, self._d_uri(*words[0:2]))
            ])

    def _siphon_parentdr(self, triples):
        for words in self._read(os.path.join("rri", "ready", "residues.txt")):
            assert len(words) == 3
            triples.extend([
                    (self._d_uri(*words[0:2]), O.YIP_PARENT_OF, self._r_uri(*words[0:3]))
            ])

    def _siphon_aliases(self, triples):
        from urllib2 import urlopen
        BOXES = (
            #(  0,  75, "GD"),           # Gene designation
            ( 75,  95, "OLN"),          # Ordered locus name
            ( 95, 106, "SP_ACC"),       # Swiss-Prot accession
            #(106, 118, "SP_NAME"),      # Swiss-Prot entry name
            #(118, 128, "SGD_ACCESSION"),# SGD accession
            # XXX there are three more fields we do not care about
        )
        try:
            url = urlopen("http://www.uniprot.org/docs/yeast.txt")
            state = 0
            for line in url.read().split("\n"):
                if "____" in line or "----" in line:
                    state += 1
                elif state == 5:
                    parts = { box[2]: map(str.strip, line[box[0]:box[1]].strip().split(";")) for box in BOXES }
                    assert len(parts["OLN"]) == 1
                    assert len(parts["SP_ACC"]) == 1
                    yip_acc = O.uri(O.YIP_PROTEIN, parts["OLN"][0])
                    sp_acc = U(O.UNIPROT_ID + parts["SP_ACC"][0])
                    triples.append((yip_acc, O.OWL.sameAs, sp_acc))
            url.close()
        except Exception, e:
            print e

    def get_entities(self):
        ps, ds, rs = [], [], []
        for words in self._read(os.path.join("ppi", "ready", "proteins.txt")):
            ps.append(words[0])
        for words in self._read(os.path.join("ddi", "ready", "domains.txt")):
            ds.append("_".join(words))
        for words in self._read(os.path.join("rri", "ready", "residues.txt")):
            rs.append("_".join(words))
        return ps, ds, rs

    def get_pairs(self):
        pps = [(words[0], words[1]) for words
               in (self._read(os.path.join("ppi", "ready", "goldPosProteinPairs.txt")) + \
                   self._read(os.path.join("ppi", "ready", "goldNegProteinPairs.txt")))]
        dds = [(words[0:2], words[2:4]) for words
               in (self._read(os.path.join("ddi", "ready", "goldPosDomainPairs.txt")) + \
                   self._read(os.path.join("ddi", "ready", "goldNegDomainPairs.txt")))]
        rrs = [(words[0:3], words[3:6]) for words
               in (self._read(os.path.join("rri", "ready", "goldPosResiduePairs.txt")) + \
                   self._read(os.path.join("rri", "ready", "goldNegResiduePairs.txt")))]
        return pps, dds, rrs
