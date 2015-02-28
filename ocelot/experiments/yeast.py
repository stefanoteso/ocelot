# -*- coding: utf-8 -*-

from .base import _Experiment

import os
import numpy as np
from ocelot.services import _cls, CDHit

class YeastExperiment(_Experiment):
    """New experiment based on SGD and iPfam.

    This experiment is structured exactly the same as the Yip et al. [Yip09]_
    experiment, but performed on a different dataset based on SGD, BioGRID
    and iPfam.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    # TODO compute composition features
    # TODO compute complexity features
    # TODO compute conservation (profile) features
    # TODO compute secondary features
    # TODO compute cell-cycle gene expression features (correlations)
    # TODO compute environment-response gene expression features (correlations)
    # TODO read in the Y2H raw data
    # TODO read in the TAP-MS raw data
    def __init__(self, *args, **kwargs):
        super(YeastExperiment, self).__init__(*args, **kwargs)

    def _query_proteins_and_info(self):
        """Queries the ORFs from SGD along with their sequences, etc."""
        ans = self.query("""
        SELECT ?orf ?seq
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_sequence ?seq .
        }}
        """)
        assert ans and len(ans) and "results" in ans
        p_to_seq = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p = bindings[u"orf"].split(".")[-1]
            p_to_seq[p] = bindings[u"seq"]

        ans = self.query("""
        SELECT ?orf ?fun
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_goslim ?fun .
        }}
        """)
        assert ans and len(ans) and "results" in ans
        p_to_fun = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p = bindings[u"orf"].split(".")[-1]
            fun = bindings[u"fun"].split("#")[-1]
            if not p in p_to_fun:
                p_to_fun[p] = set([ fun ])
            else:
                p_to_fun[p].update([ fun ])

        p_to_feat = {}
        return set(p_to_seq.keys()), p_to_seq, p_to_fun, p_to_feat

    def _get_proteins_and_info(self):
        """Tries to load the SGD ORFs from disk, or falls back to SPARQL."""
        # XXX not exactly elegant
        try:
            assert not self.force_update, "forcing update"
            ps          = self._depickle("ps.txt")
            p_to_seq    = self._depickle("p_to_seq.txt")
            p_to_fun    = self._depickle("p_to_fun.txt")
            p_to_feat   = self._depickle("p_to_feat.txt")
        except Exception, e:
            print _cls(self), ":", e
            ps, p_to_seq, p_to_fun, p_to_feat = \
                self._query_proteins_and_info()
            self._pickle(ps, "ps.txt")
            self._pickle(p_to_seq, "p_to_seq.txt")
            self._pickle(p_to_fun, "p_to_fun.txt")
            self._pickle(p_to_feat, "p_to_feat.txt")
        return ps, p_to_seq, p_to_fun, p_to_feat

    def _get_protein_interactions(self, manual_only = False):
        query = """
        SELECT ?orf1 ?orf2
        FROM <{default_graph}>
        WHERE {{
            ?orf1 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                owl:sameAs ?feat1 .
            ?orf2 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                owl:sameAs ?feat2 .
            ?feat1 a ocelot:sgd_feature .
            ?feat2 a ocelot:sgd_feature .
            ?int a ocelot:sgd_int ;
                ocelot:sgd_int_has_bait ?feat1 ;
                ocelot:sgd_int_has_hit ?feat2 ;
        """
        if manual_only:
            query += "      ocelot:sgd_int_has_cur_type ocelot:sgd_int_cur_type.manually_curated ;"
        query += """
                ocelot:sgd_int_has_type ocelot:sgd_int_type.physical_interactions .
        }}"""
        ans = self.query(query)
        assert ans and len(ans) and "results" in ans
        pp_pos = set()
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p1 = bindings[u"orf1"].split(".")[-1]
            p2 = bindings[u"orf2"].split(".")[-1]
            pp_pos.update([ (p1,p2) ])
        return pp_pos

    def _get_negative_protein_interactions(self, ps, pos):
        """Samples the negative interactions from the complement graph.

        Uses the sampling method described in [Yu10]_.
        """
        neighbors_of = { p: set() for p in ps }
        for p, q in pos:
            neighbors_of[p].update([q])
            neighbors_of[q].update([p])
        for p in neighbors_of:
            neighbors_of[p] = list(neighbors_of[p])
        ps_degree = sorted([(len(neighbors_of[p]), p) for p in ps ],
                            reverse = True)

        def is_candidate(p, q, neg):
            # TODO subtract BioGRID (or even better STRING)
            return p != q and not (p, q) in pos and not (p, q) in neg

        neg = set()
        for i, (degree, p) in enumerate(ps_degree):
            # Figure out all candidate negative interactions for ``p``, then
            # pick exactly ``degree`` of them uniformly at random; if not
            # enough negatives are available, add as many as possible.
            candidate_neg = filter(lambda pq: is_candidate(pq[0], pq[1], neg),
                                [(p, q) for q in ps if not q in neighbors_of[p]])
            num_samples = min(degree, len(candidate_neg))
            if num_samples == 0:
                # Proteins are sorted by degree, so we can break here
                break
            print _cls(self), "| {}/{}, '{}: sampling {}/{} negatives (actually {})" \
                .format(i+1, len(ps), p, degree, len(candidate_neg), num_samples)
            if num_samples < degree:
                print _cls(self), "| Warning: not enough negative candidates!"
            sample = [candidate_neg[i] for i
                      in np.random.permutation(len(candidate_neg))[:num_samples]]
            # Make sure that the negative interactions are symmetrical
            neg.update(sample)
            neg.update([(q, p) for p, q in sample])
        return neg

    def _get_folds(ps, pos, neg):
        # First, the known protein interactions in iPfam were divided into ten
        # folds. Then, each domain instance interaction and each residue
        # interaction was put into the fold in which the parent protein
        # interaction was assigned.  Finally, the remaining protein
        # interactions and all the negative sets were randomly divided into ten
        # folds.
        pass

    @staticmethod
    def _check_p_pp_are_sane(ps, pps):
        assert all((p, q) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"

    def run(self):
        """Run the yeast prediction experiment."""
        ps, p_to_seq, p_to_fun, p_to_feat = self._get_proteins_and_info()
        print _cls(self), ": found {} proteins".format(len(ps))

        # Filter out sequences shorter than 30 residues
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= 30, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), 30)

        # Cluster sequences with CD-HIT
        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, filtered_p_clusters = CDHit().run(filtered_p_seq, threshold = 0.8)
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(filtered_p_clusters), 0.8)

        # Query the hq protein-protein interactions
        pp_pos_hq = self._get_protein_interactions(manual_only = True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality protein-protein interactions (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_p_pp_are_sane(ps, pp_pos_hq)

        # Query the hq+lq protein-protein interactions
        pp_pos_lq = self._get_protein_interactions(manual_only = False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality protein-protein interactions (density = {})" \
                .format(len(pp_pos_lq), density)
        self._check_p_pp_are_sane(ps, pp_pos_lq)

        # Here we pass the low-quality interactions so as to get a better
        # approximation of the negative set.
        pp_neg = self._get_negative_protein_interactions(ps, pp_pos_lq)
        print _cls(self), ": sampled {} pp negative interactions" \
                .format(len(pp_neg))

        # TODO retrieve dom-dom and res-res interaction instances

        # TODO create the folds, with redundancy reduced training sets
        pp_folds = self._get_folds(ps, pp_pos_hq, pp_neg)
