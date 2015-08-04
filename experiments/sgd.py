# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools as it
from collections import defaultdict
from ocelot.services import _cls, CDHit
from ocelot.go import GODag, GOTerm
from ocelot.kernels import *
from . import _Experiment, Stage
from .yeast import *

class SGDExperiment(_Experiment):
    """New experiment based on SGD and iPfam.

    This experiment is structured exactly the same as the Yip et al. [Yip09]_
    experiment, but performed on a different dataset based on SGD, BioGRID
    and iPfam.

    .. note::

        Somewhat unrelated. For up-to-date statistics on GO annotations,
        see `<http://geneontology.org/page/current-go-statistics>`_.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    def __init__(self, *args, **kwargs):
        self._go_aspects = kwargs.get("go_aspects", None)
        self._max_go_depth = kwargs.get("max_go_depth", None)
        self._min_go_annot = kwargs.get("min_go_annot", None)
        super(SGDExperiment, self).__init__(*args, **kwargs)

    def _get_sgd_ids(self):
        query = """
        SELECT ?orf ?seq
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified .
        }}
        """
        ids = set(bindings[u"orf"].split(".")[-1]
                  for bindings in self.iterquery(query, n = 1))
        return sorted(list(ids))

    def _get_sgd_id_to_seq(self):
        query = """
        SELECT ?orf ?seq
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
                 ocelot:sgd_id_has_sequence ?seq .
        }}
        """
        sgd_id_to_seq = {}
        for bindings in self.iterquery(query, n = 2):
            sgd_id = bindings[u"orf"].split(".")[-1]
            sgd_id_to_seq[sgd_id] = bindings[u"seq"]
        return sgd_id_to_seq

    def _get_sgd_id_to_term_ids(self):
        query = """
        SELECT ?orf ?fun
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
                 ocelot:sgd_id_has_goslim ?fun .
        }}
        """
        sgd_id_to_fun = {}
        for bindings in self.iterquery(query, n = 2):
            sgd_id = bindings[u"orf"].split(".")[-1]
            fun = bindings[u"fun"].split("#")[-1]
            if not sgd_id in sgd_id_to_fun:
                sgd_id_to_fun[sgd_id] = set([ fun ])
            else:
                sgd_id_to_fun[sgd_id].add(fun)
        return sgd_id_to_fun

    def _get_sgd_id_to_feat(self):
        query = """
        SELECT ?orf ?feat
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
                 owl:sameAs ?feat .
            ?feat a ocelot:sgd_feature .
        }}
        """
        sgd_id_to_feat = {}
        for bindings in self.iterquery(query, n = 2):
            sgd_id = bindings[u"orf"].split(".")[-1]
            feat   = bindings[u"feat"].split(".")[-1]
            if not sgd_id in sgd_id_to_feat:
                sgd_id_to_feat[sgd_id] = set([ feat ])
            else:
                sgd_id_to_feat[sgd_id].add(feat)
        for sgd_id, featset in sgd_id_to_feat.iteritems():
            assert len(featset) == 1
            sgd_id_to_feat[sgd_id] = list(featset)[0]
        return sgd_id_to_feat

    def _get_sgd_id_to_context(self):
        query = """
        SELECT ?orf ?chrom ?strand ?start ?stop
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified .
            ?orf ocelot:sgd_id_in_chromosome ?chrom .
            ?orf ocelot:sgd_id_in_strand ?strand .
            ?orf ocelot:sgd_id_starts_at ?start .
            ?orf ocelot:sgd_id_stops_at ?stop .
        }}
        """
        p_to_context = {}
        for bindings in self.iterquery(query, n = 5):
            orf     = bindings[u"orf"].split(".")[-1]
            chrom   = bindings[u"chrom"].split(".")[-1]
            strand  = bindings[u"strand"]
            start   = bindings[u"start"]
            stop    = bindings[u"stop"]
            assert strand in ("C", "W")
            p_to_context[orf] = \
                (chrom, min(start, stop) + 0.5 * np.fabs(start - stop))
        return p_to_context

    def _get_sgd_pin(self, ps = None, manual_only = False):
        query = """
        SELECT DISTINCT ?orf1 ?orf2
        FROM <{default_graph}>
        WHERE {{
            ?orf1 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
                owl:sameAs ?feat1 .
            ?orf2 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
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
        pp_pos = set()
        for bindings in self.iterquery(query, n = 2):
            p1 = bindings[u"orf1"].split(".")[-1]
            p2 = bindings[u"orf2"].split(".")[-1]
            pp_pos.update([(p1,p2), (p2,p1)])
        if not ps is None:
            # Filter out all pairs that do not fall in ``ps``
            ps, filtered_pp_pos = set(ps), set()
            for (p1, p2) in pp_pos:
                if p1 in ps and p2 in ps:
                    filtered_pp_pos.add((p1, p2))
            pp_pos = filtered_pp_pos
        return pp_pos

    def _get_string_pin(self, ps = None):
        query = """
        SELECT DISTINCT ?orf1 ?orf2
        FROM <{default_graph}>
        WHERE {{
            ?orf1 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified .
            ?orf2 a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified .
            ?id1 a ocelot:STRING_ID ;
                owl:sameAs ?orf1 .
            ?id2 a ocelot:STRING_ID ;
                owl:sameAs ?orf2 .
            ?id1 ?mode ?id2 .
        }}"""
        pp_pos = set()
        for bindings in self.iterquery(query, n = 2):
            p1 = bindings[u"orf1"].split(".")[-1]
            p2 = bindings[u"orf2"].split(".")[-1]
            pp_pos.update([(p1,p2), (p2,p1)])
        if not ps is None:
            # Filter out all pairs that do not fall in ``ps``
            ps, filtered_pp_pos = set(ps), set()
            for (p1, p2) in pp_pos:
                if p1 in ps and p2 in ps:
                    filtered_pp_pos.add((p1, p2))
            pp_pos = filtered_pp_pos
        return pp_pos

    @staticmethod
    def _get_neighbors_and_degree(ps, pps):
        neighbors_of = { p: set() for p in ps }
        for p, q in pps:
            neighbors_of[p].add(q)
            neighbors_of[q].add(p)
        return neighbors_of

    @staticmethod
    def _check_p_pp_are_sane(ps, pps):
        """Checks that the protein pairs are (i) symmetric, and (ii) entirely
        contained in the list of proteins."""
        assert all((q, p) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"

    def _get_sgd_pins(self, ps):
        """Computes the high-quality positive interactions, high+low-quality
        positive interactions, and the negative interactions."""
        pp_pos_hq = self._get_sgd_pin(ps=ps, manual_only=True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_p_pp_are_sane(ps, pp_pos_hq)

        # Query all (high+low-quality) protein-protein interactions
        pp_pos_lq = self._get_sgd_pin(ps=ps, manual_only=False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality PPIs (density = {})" \
                .format(len(pp_pos_lq), density)
        self._check_p_pp_are_sane(ps, pp_pos_lq)

        # Query all (literally) protein-protein actions annotated in STRING
        pp_pos_string = self._get_string_pin(ps=ps)
        density = float(len(pp_pos_string)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} STRING-quality PPIs (density = {})" \
                .format(len(pp_pos_string), density)
        self._check_p_pp_are_sane(ps, pp_pos_string)

        return pp_pos_hq, pp_pos_lq | pp_pos_string

    @staticmethod
    def _desymmetrize(pairs):
        pairs_asymm = set()
        for p, q in pairs:
            if not (q, p) in pairs:
                pairs_asymm.add((p, q))
        return pairs_asymm

    @staticmethod
    def _symmetrize(pairs):
        return set(pairs) | set((q, p) for p, q in pairs)

    def _compute_negative_pin(self, ps, pp_pos_hq, pp_pos_lq):

        # Take the complement of the symmetrical positives, subtract the known
        # positives (both high-quality and low-quality, just to be sure); the
        # complement will be symmetrical, since the HQ and LQ interactions are
        print _cls(self), ": computing the complement of the positive PIN"
        pp_neg_all = set(it.product(ps, ps)) - (pp_pos_lq | pp_pos_hq)

        print _cls(self), ": computing asymmetric negative PIN"
        pp_neg_asymm = list(self._desymmetrize(pp_neg_all))

        print _cls(self), ": computing asymmetric positive PIN"
        pp_pos_asymm = self._desymmetrize(pp_pos_hq)

        del pp_pos_hq
        del pp_pos_lq
        del pp_neg_all

        # Sample the negative interactions from the complement
        print _cls(self), ": sampling the negative PIN"
        pi = self._rng.permutation(len(pp_neg_asymm))
        pp_neg_half = set(pp_neg_asymm[pi[i]] for i in xrange(len(pp_pos_asymm)))

        # Symmetrize the negatives
        pp_neg = self._symmetrize(pp_neg_half)
        self._check_p_pp_are_sane(ps, pp_neg)

        return pp_neg,

    def _get_sgd_din(self):
        query = """
        SELECT DISTINCT ?id ?family ?chain ?evalue ?complex
        FROM <{default_graph}>
        WHERE {{
            # Pick all SGD ORFs
            ?orf a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified ;
                owl:sameAs ?feat .

            # Map them to SGD feature names
            ?feat a ocelot:sgd_feature .

            # Map them to PDB chains via SGD's precomputed homology mappings
            ?_h a ocelot:sgd_pdb_homology ;
                ocelot:sgd_pdb_has_query ?feat ;
                ocelot:sgd_pdb_has_target ?chain ;
                ocelot:sgd_pdb_alignment ?evalue .

            # Find all iPfam domain instances including that chain
            ?region a ocelot:ipfam_region ;
                ocelot:ipfam_region_instance_of ?family ;
                ocelot:ipfam_region_occurs_in ?chain .

            ?_ri a ocelot:ipfam_region_int ;
                ocelot:ipfam_region_int_has_region ?region ;
                ocelot:ipfam_region_int_occurs_in ?complex .

            # TODO filter out all non-yeast complexes
        }}"""
        dd_pos = set()
        for bindings in self.iterquery(query, n=5):
            print bindings
        sys.exit(1)
        return dd_pos

    def _compute_p_colocalization_kernel(self, ps):
        p_to_context = self._get_sgd_id_to_context()
        contexts = [p_to_context[p] for p in ps]
        return self._compute_kernel(ColocalizationKernel, contexts, gamma=1.0)

    def _compute_p_gene_expression_kernel(self, ps, p_to_feat):
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}
        return self._compute_kernel(SGDGeneExpressionKernel,
                                    feat_to_i, self.src, 
                                    ["Gasch_2000_PMID_11102521",
                                     "Spellman_1998_PMID_9843569"])

    def _compute_p_complex_kernel(self, ps, p_to_feat):
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}
        return self._compute_kernel(YeastProteinComplexKernel,
                                    feat_to_i, self.src)

    def _compute_p_interpro_kernel(self, ps):
        return self._compute_kernel(InterProKernel, ps, self.dst,
                                    use_evalue=False)

    def _compute_p_pssm_kernel(self, ps, p_to_seq):
        return self._compute_kernel(PSSMKernel, ps, p_to_seq, self.dst,
                                    k=4, threshold=6.0)

    def _compute_pp_kernel(self, ps, folds, submatrix):
        p_to_i = {p: i for i, p in enumerate(ps)}

        pps = set()
        for fold in folds:
            pps.update((p1, p2) for p1, p2, _ in fold)
        pps = sorted(pps)

        pp_indices = [(p_to_i[p1], p_to_i[p2]) for p1, p2 in pps]
        return self._compute_kernel(PairwiseKernel, pp_indices, submatrix)

    def _permute(self, l):
        pi = list(self._rng.permutation(len(l)))
        return [l[pi[i]] for i in xrange(len(l))]

    @staticmethod
    def _check_folds_are_sane(folds):

        # Check that interactions are symmetric
        for k, fold in enumerate(folds):
            assert all((q, p, state) in fold for (p, q, state) in fold), "fold {} is not symmetric".format(k)

        # Check that folds do not overlap
        for (k1, fold1), (k2, fold2) in it.product(enumerate(folds), enumerate(folds)):
            if k1 >= k2:
                continue
            for p1, p2, state1 in fold1:
                for q1, q2, state2 in fold2:
                    assert (p1, p2) != (q1, q2), "folds {} and {} overlap ({} overlaps with {})".format(k1, k2, (p1, p2, state1), (q1, q2, state2))

        # Check that no pair occurs in both states
        interactions = set()
        for fold in folds:
            interactions.update(fold)
        for p1, p2, state in interactions:
            assert not (p1, p2, not state) in interactions, "folds are inconsistent"

    def _print_fold_quality(self, folds, pp_to_terms, all_pp_terms, term_to_i, average):
        # Evaluate the fold quality
        for k, fold in enumerate(folds):
            counts = np.zeros((len(all_pp_terms),))
            for p1, p2, _ in fold:
                for term in pp_to_terms[(p1, p2)]:
                    counts[term_to_i[term]] += 1
            print _cls(self), " fold{}: {} pairs, L1 distance from average term counts = {}" \
                .format(k, len(fold), np.linalg.norm((average - counts) / len(all_pp_terms), ord = 1))

    def _compute_folds(self, pp_pos_hq, pp_pos_lq, p_to_terms, num_folds=10):
        """Generates the folds.

        The interactions are split into folds according to their functions.
        No two protein-protein interactions can occur in two distinct folds.

        :param pp_pos_hq: high-quality positive interactions.
        :param pp_pos_lq: low-quality positive interactions, if any.
        :param p_to_terms: map between proteins and ``GOTerm``'s.
        :param num_folds: number of folds to generate (default: ``10``).
        :returns: a list of sets of protein pairs.
        """
        assert isinstance(pp_pos_hq, set)
        assert isinstance(pp_pos_lq, set)

        # Map from protein *pairs* to the GO terms that either protein is
        # annotated with, and vice-versa
        pp_to_terms = {(p1,p2): (set(p_to_terms[p1]) | set(p_to_terms[p2]))
                       for p1, p2 in pp_pos_hq}

        term_to_pps = defaultdict(set)
        for pp, terms in pp_to_terms.iteritems():
            for term in terms:
                term_to_pps[term].add(pp)

        # Gather all GO terms and assign them an index
        all_pp_terms = term_to_pps.keys()
        term_to_i = {term: i for i, term in enumerate(all_pp_terms)}

        # Compute the ideal (perfect average) GO term distribution
        average = np.zeros((len(all_pp_terms),))
        for pp, terms in pp_to_terms.iteritems():
            for term in terms:
                average[term_to_i[term]] += 1
        average *= 1.0 / num_folds

        # Sanity check. Note that self-interactions are fine.
        for term, term_pps in term_to_pps.iteritems():
            for p1, p2 in term_pps:
                assert (p2, p1) in term_pps

        # Generate the folds
        print _cls(self), ": generating the folds..."

        folds = [set() for _ in range(num_folds)]
        while len(term_to_pps):
            # Pick the term with the least unprocessed protein-protein pairs
            cur_term, cur_pps = min(term_to_pps.iteritems(),
                                    key = lambda term_pps: len(term_pps[1]))
            print _cls(self), ": best term: {} (num pairs = {})".format(cur_term, len(cur_pps))

            # Evenly distribute the associated protein-protein pairs among all
            # folds, taking into account the fact that if (p1, p2) is in
            # cur_pps, then (p2, p1) is in cur_pps as well. Also, make sure
            # to allow self-interactions.
            pps_to_add = set()
            for p1, p2 in cur_pps:
                if not (p2, p1) in pps_to_add:
                    pps_to_add.add((p1, p2))
            folds = self._permute(folds)
            for i, (p1, p2) in enumerate(pps_to_add):
                folds[i % num_folds].update([(p1, p2, True), (p2, p1, True)])

            # Update term_to_pps
            new_term_to_pps = {}
            for term in term_to_pps:
                # Skip the term we just processed
                if term == cur_term:
                    continue
                # Skip the protein pairs annotated with it
                new_term_pps = set(pp for pp in term_to_pps[term]
                                   if not pp in cur_pps)
                if len(new_term_pps) == 0:
                    continue
                new_term_to_pps[term] = new_term_pps
            term_to_pps = new_term_to_pps

        print _cls(self), ": checking fold sanity..."
        self._check_folds_are_sane(folds)

        print _cls(self), ": fold quality (positives only):"
        self._print_fold_quality(folds, pp_to_terms, all_pp_terms, term_to_i, average)

        # Now that we partitioned all interactions into folds, let's add the
        # negatives interactions -- by sampling them at random
        print _cls(self), ": adding negative interactions..."

        candidate_pos_pps = {(p1, p2) for p1, p2 in pp_pos_hq | pp_pos_lq} | \
                            {(p2, p1) for p1, p2 in pp_pos_hq | pp_pos_lq}

        # All candadate negative pairs sampled so far, so that we do not sample
        # the same negative twice
        all_candidate_neg_pps = set()

        for fold in folds:
            ps_in_fold = set(p for p, _, _ in fold) | set(p for _, p, _ in fold)

            # Compute the candidate negative pairs by subtracting the candidate
            # positive pairs from the complement of the pairs in the fold
            candidate_neg_pps = set(it.product(ps_in_fold, ps_in_fold)) - candidate_pos_pps - all_candidate_neg_pps

            # De-symmetrize the candidate negative pairs
            temp = set()
            for p1, p2 in candidate_neg_pps:
                if not (p2, p1) in temp:
                    temp.add((p1, p2))
            candidate_neg_pps = list(temp)

            # De-symmetrize the positive fold
            temp = set()
            for p1, p2, _ in fold:
                if not (p2, p1) in temp:
                    temp.add((p1, p2))

            # Randomly sample a number of negative pairs equal to the number of
            # positives pairs (actually only half of that; we symmetrize the
            # negatives later)
            pi = self._rng.permutation(len(candidate_neg_pps))
            sampled_neg_pps = set(candidate_neg_pps[pi[i]] for i in xrange(len(temp)))

            # Update with the sampled negative pairs
            all_candidate_neg_pps.update((p1, p2) for p1, p2 in sampled_neg_pps)
            all_candidate_neg_pps.update((p2, p1) for p1, p2 in sampled_neg_pps)

            # Assemble the fold
            fold.update((p1, p2, False) for p1, p2 in sampled_neg_pps)
            fold.update((p2, p1, False) for p1, p2 in sampled_neg_pps)

        print _cls(self), ": checking fold sanity..."
        self._check_folds_are_sane(folds)

        # XXX pp_to_terms is not updated at this point, do not bother printing
        # the fold quality...

        return folds,

    def _dump_dataset_statistics(self, ps, p_to_seq, p_to_term_ids, dag, prefix):
        """Dump a few dataset statistics."""

        # Draw the histogram of protein lengths
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = [len(p_to_seq[p]) for p in ps]
        ax.hist(x, range(min(x), max(x) + 50, 50), color = "red", alpha = 0.8)
        ax.yaxis.grid(True)
        fig.savefig(os.path.join(self.dst, "p-{}-length-hist.png".format(prefix)),
                    bbox_inches = "tight", pad_inches = 0)

        # Draw the GO annotation depth
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = []
        for term in dag._id_to_term.itervalues():
            x.extend([term.level] * len(term.proteins))
        ax.hist(x, bins = 15, color = "red", alpha = 0.8)
        ax.yaxis.grid(True)
        fig.savefig(os.path.join(self.dst, "p-{}-annot-depth.png".format(prefix)),
                    bbox_inches = "tight", pad_inches = 0)

        # Print the protein->function map
        with open(os.path.join(self.dst, "p-{}-p-to-term-ids.txt".format(prefix)), "wb") as fp:
            for p, term_ids in p_to_term_ids.iteritems():
                fp.write("{}: {}\n".format(p, ",".join(sorted(term_ids))))

        # Print the function->protein map
        with open(os.path.join(self.dst, "p-{}-term-id-to-ps.txt".format(prefix)), "wb") as fp:
            for term_id, term in dag._id_to_term.iteritems():
                fp.write("{}: {}\n".format(term_id, ",".join(sorted(term.proteins))))

        annotated_term_ids = [term_id for term_id, term in dag._id_to_term.iteritems()
                              if len(term.proteins)]
        annotated_term_id_to_i = {term_id: i for i, term_id
                                  in enumerate(annotated_term_ids)}

        # Plot function co-occurrence
        tt_condprob = []
        for term_id1, term_id2 in it.product(annotated_term_ids, annotated_term_ids):
            ps1 = dag._id_to_term[term_id1].proteins
            ps2 = dag._id_to_term[term_id2].proteins
            condprob = len(ps1 & ps2) / float(len(ps1))
            tt_condprob.append((term_id1, term_id2, condprob))
        tt_condprob = sorted(tt_condprob, key = lambda triple: triple[-1],
                               reverse = True)
        with open(os.path.join(self.dst, "p-{}-term-condprob.txt".format(prefix)), "wb") as fp:
            for term_id1, term_id2, condprob in tt_condprob:
                fp.write("{},{}: {}\n".format(term_id1, term_id2, condprob))

        matrix = np.zeros((len(annotated_term_ids), len(annotated_term_ids)))
        for term_id1, term_id2, condprob in tt_condprob:
            i = annotated_term_id_to_i[term_id1]
            j = annotated_term_id_to_i[term_id2]
            matrix[i,j] = condprob

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.matshow(matrix, interpolation = "nearest", cmap = cm.OrRd)
        fig.savefig(os.path.join(self.dst, "p-{}-term-condprob.png".format(prefix)),
                    bbox_inches="tight", pad_inches=0)

        # Print consistency w.r.t. mutual exclusion of children
        tt_mutex = []
        for term_id in annotated_term_ids:
            term = dag._id_to_term[term_id]
            for (child_id1, rel1), (child_id2, rel2) in it.product(term.children, term.children):
                if child_id1 == child_id2:
                    continue
                if rel1 != "is_a" or rel2 != "is_a":
                    continue
                ps1 = dag._id_to_term[child_id1].proteins
                ps2 = dag._id_to_term[child_id2].proteins
                ni, nu = len(ps1 & ps2), len(ps1 | ps2)
                tt_mutex.append((child_id1, child_id2, nu, ni))
        tt_mutex = sorted(tt_mutex, key = lambda quad: quad[-1], reverse = True)
        with open(os.path.join(self.dst, "p-{}-term-mutex.txt".format(prefix)), "wb") as fp:
            for term_id1, term_id1, nu, ni in tt_mutex:
                fp.write("{},{}: {}/{}={}\n".format(term_id1, term_id2, ni, nu, ni / float(nu) if nu > 0 else -1.0))

        # TODO: plot GO consistency w.r.t. interactions

    def _load_sgd_resources(self):
        """Retrieves all SGD-derived resources.

        Namely:

        - the list of protein (validated ORF) IDs.
        - the ID -> yeast protein name map.
        - the ID -> sequence map.
        - the ID -> GO term annotations map.

        :returns: the four items above.
        """
        ps              = self._get_sgd_ids()
        p_to_feat       = self._get_sgd_id_to_feat()
        p_to_seq        = self._get_sgd_id_to_seq()
        p_to_term_ids   = self._get_sgd_id_to_term_ids()

        print _cls(self), ": found {} proteins".format(len(ps))
        for p in ps:
            assert p in p_to_seq, "'{}' has no sequence".format(p)
            assert p in p_to_term_ids, "'{}' has no GO annotation".format(p)
            assert p in p_to_feat, "'{}' has no associated feature ID".format(p)

        return ps, p_to_feat, p_to_seq, p_to_term_ids

    def _filter_ps(self, ps, p_to_seq, min_sequence_len, cdhit_threshold):
        """Filter out short sequences and cluster them with CD-HIT."""

        # Filter out sequences shorter than min_sequence_len residues
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= min_sequence_len, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), min_sequence_len)

        # Cluster the filtered sequences with CD-HIT
        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, clusters = CDHit().run(filtered_p_seq, threshold = cdhit_threshold)
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(clusters), cdhit_threshold)

        # Take one protein from each cluster
        filtered_ps = [list(cluster)[0][0] for cluster in clusters]
        print _cls(self), ": there are {} filtered proteins".format(len(filtered_ps))

        return filtered_ps,

    def _load_go_and_filter(self, filtered_ps, p_to_term_ids):
        """Load the GO OBO file, fill it in, and prune it."""

        # Load the GO OBO file
        dag = GODag(os.path.join(self.src, "GO", "go-basic.obo"))

        # Fill in the GO data structure with the protein annotations and
        # propagate them to the root
        propagated_p_to_term_ids = dag.annotate(p_to_term_ids, propagate=True)

        print _cls(self), ": '{}' annotations in the DAG" \
                            .format(sum(len(dag._id_to_term[term_id].proteins)
                                        for term_id in dag._id_to_term.iterkeys()))
        tot_increase = 0.0
        for p in p_to_term_ids:
            tot_increase += len(propagated_p_to_term_ids[p]) / float(len(p_to_term_ids[p]))

        print _cls(self), ": propagated GO annotations, average increase is '{}'" \
                            .format(tot_increase / len(p_to_term_ids))

        aspect_to_ps = dag.get_proteins_by_aspect()
        print _cls(self), ": annotations by aspect: {}".format([(aspect, len(ps)) for aspect, ps in aspect_to_ps.iteritems()])

        filtered_p_to_term_ids = dag.preprocess(filtered_ps,
                                                aspects = self._go_aspects,
                                                min_annot = self._min_go_annot,
                                                max_depth = self._max_go_depth)

        return dag, filtered_p_to_term_ids

    def run(self, min_sequence_len=50, cdhit_threshold=0.75, dump_stats=False):
        """Run the yeast prediction experiment.

        :param min_sequence_len: minimum sequence length (default: ``50``).
        :param cdhit_threshold: CD-HIT threshold (default: ``0.75``).
        :param dump_stats: whether to dump statistics (default: ``False``).
        """
        STAGES = (

            Stage(self._load_sgd_resources,
                  [], ['ps', 'p_to_feat', 'p_to_seq', 'p_to_term_ids']),

            Stage(self._filter_ps,
                  ['ps', 'p_to_seq', 'min_sequence_len', 'cdhit_threshold'],
                  ['filtered_ps']),

            Stage(self._compute_p_colocalization_kernel,
                  ['filtered_ps'], ['p_colocalization_kernel']),

            Stage(self._compute_p_gene_expression_kernel,
                  ['filtered_ps', 'p_to_feat'], ['p_gene_expression_kernel']),

            Stage(self._compute_p_complex_kernel,
                  ['filtered_ps', 'p_to_feat'], ['p_complex_kernel']),

            Stage(self._compute_p_interpro_kernel,
                  ['filtered_ps'], ['p_interpro_kernel']),

            Stage(self._compute_p_pssm_kernel,
                  ['filtered_ps', 'p_to_seq'], ['p_pssm_kernel']),

            Stage(self._load_go_and_filter,
                  ['filtered_ps', 'p_to_term_ids'],
                  ['filtered_dag', 'filtered_p_to_term_ids']),

            Stage(self._get_sgd_pins,
                  ['filtered_ps'], ['pp_pos_hq', 'pp_pos_lq']),

            Stage(self._compute_negative_pin,
                  ['filtered_ps', 'pp_pos_hq', 'pp_pos_lq'],
                  ['pp_neg']),

            Stage(self._compute_folds,
                  ['pp_pos_hq', 'pp_pos_lq', 'filtered_p_to_term_ids'], ['folds']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'folds', 'p_colocalization_kernel'],
                  ['pp_colocalization_kernel']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'folds', 'p_gene_expression_kernel'],
                  ['pp_gene_expression_kernel']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'folds', 'p_complex_kernel'],
                  ['pp_complex_kernel']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'folds', 'p_interpro_kernel'],
                  ['pp_interpro_kernel']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'folds', 'p_pssm_kernel'],
                  ['pp_pssm_kernel']),

            Stage(lambda *args, **kwargs: None,
                  ['pp_colocalization_kernel'],
                  ['pp_kernels']),
        )

        TARGETS = (
            'pp_neg',
        )

        context = {
            "min_sequence_len": min_sequence_len,
            "cdhit_threshold": cdhit_threshold,
            "dump_states": dump_stats,
        }
        self._make(STAGES, TARGETS, context=context)
