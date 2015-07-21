# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools as it
from ocelot.services import _cls, CDHit
from ocelot.go import GODag, GOTerm
from ocelot.kernels import *
from . import _Experiment
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

    def _get_sgd_id_to_fun(self):
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

    def _get_negative_pin(self, ps, pos):
        """Samples the negative interactions from the complement graph.

        Uses the sampling method described in [Yu10]_.
        """
        def is_candidate(p, q, neg):
            # TODO subtract BioGRID (or even better STRING)
            return p != q and not (p, q) in pos and not (p, q) in neg

        neighbors_of = self._get_neighbors_and_degree(ps, pos)
        ps_degree = sorted([(len(neighbors_of[p]), p) for p in ps ],
                            reverse = True)

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
        pp_pos_hq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_hq",
                                 ps = ps, manual_only = True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_p_pp_are_sane(ps, pp_pos_hq)

        # Query all (high+low-quality) protein-protein interactions
        pp_pos_lq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_lq",
                                 ps = ps, manual_only = False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality PPIs (density = {})" \
                .format(len(pp_pos_lq), density)
        self._check_p_pp_are_sane(ps, pp_pos_lq)

        # Query all (literally) protein-protein actions annotated in STRING
        pp_pos_string = self._cached(self._get_string_pin,
                                     "sgd_id_interactions_string",
                                     ps = ps)
        density = float(len(pp_pos_string)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} STRING-quality PPIs (density = {})" \
                .format(len(pp_pos_string), density)
        self._check_p_pp_are_sane(ps, pp_pos_string)

        return pp_pos_hq, pp_pos_lq | pp_pos_string

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

    def _compute_protein_kernels(self, ps, p_to_feat, p_to_seq, folds):
        """Computes all the kernels and pairwise kernels for proteins.

        The order in which objects are passed in is preserved.

        :param ps: list of protein ORF IDs.
        :param p_to_feat: map from protein ORF ID to protein feature ID.
        :param p_to_seq: map from protein ORF ID to protein sequence.
        :param folds: list of folds, each a list (p1, p2, label) pairs.
        """
        pps = set()
        for fold in folds:
            pps.update((p1, p2) for p1, p2, _ in fold)
        pps = sorted(pps)

        p_to_i = {p: i for i, p in enumerate(ps)}
        pp_indices = [(p_to_i[p1], p_to_i[p2]) for p1, p2 in pps]
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}

        p_to_context = self._cached(self._get_sgd_id_to_context,
                                    "sgd_id_to_context")
        contexts = [p_to_context[p] for p in ps]

        self._cached_kernel(ColocalizationKernel, len(ps),
                            "p-colocalization-kernel",
                            contexts, gamma = 1.0,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(SGDGeneExpressionKernel, len(ps),
                            "p-gene-expression-kernel",
                            feat_to_i, self.src,
                            ["Gasch_2000_PMID_11102521",
                             "Spellman_1998_PMID_9843569"],
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(YeastProteinComplexKernel, len(ps),
                            "p-complex-kernel",
                            feat_to_i, self.src,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-kernel", ps, self.dst,
                            use_evalue = False,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-score-kernel", ps, self.dst,
                            use_evalue = True,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(PSSMKernel, len(ps),
                            "p-profile-kernel", ps, p_to_seq, self.dst,
                            k = 4, threshold = 6.0,
                            do_pairwise = True, pairs = pp_indices)

    @staticmethod
    def _permute(l):
        pi = list(np.random.permutation(len(l)))
        return [l[pi[i]] for i in xrange(len(l))]

    def _get_folds(self, pp_pos, p_to_terms, num_folds = 10):
        """Generates the folds.

        The interactions are split into folds according to their functions.
        No two protein-protein interactions can occur in two distinct folds.

        :param pp_pos: list of interactions.
        :param p_to_terms: map between proteins and ``GOTerm``'s.
        :param num_folds: number of folds to generate (default: ``10``).
        :returns: a list of sets of protein pairs.
        """
        # Map from protein pairs to the GO terms that either protein is
        # annotated with
        pp_to_terms = {(p1,p2): (set(p_to_terms[p1]) | set(p_to_terms[p2]))
                       for p1, p2 in pp_pos}

        # Map from GO terms to protein pairs
        term_to_pps = {}
        for (p1, p2), terms in pp_to_terms.iteritems():
            for term in terms:
                if not term in term_to_pps:
                    term_to_pps[term] = set([(p1, p2)])
                else:
                    term_to_pps[term].add((p1, p2))

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
        folds = [set() for _ in range(num_folds)]
        while len(term_to_pps):
            # Pick the term with the least unprocessed protein-protein pairs
            cur_term, cur_pps = min(term_to_pps.iteritems(),
                                    key = lambda term_pps: len(term_pps[1]))
            print "best term: {} (num pairs = {})".format(cur_term, len(cur_pps))
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
                folds[i % num_folds].update([(p1, p2), (p2, p1)])

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

        # Sanity checks
        for fold in folds:
            for p1, p2 in fold:
                assert (p2, p1) in fold
        for (k1, fold1), (k2, fold2) in it.product(enumerate(folds), enumerate(folds)):
            assert k1 == k2 or len(fold1 & fold2) == 0

        # Evaluate the fold quality
        for k, fold in enumerate(folds):
            counts = np.zeros((len(all_pp_terms),))
            for pp in fold:
                for term in pp_to_terms[pp]:
                    counts[term_to_i[term]] += 1
            print " |fold{}| = {}, l1 distance from perfect GO term counts = {}" \
                .format(k, len(fold), np.linalg.norm((average - counts), ord = 1))

        return folds

    def _add_negatives(self, pos_folds, pp_not_neg):
        """Computes folds where negatives are randomly sampled from the
        non-positives.

        :param pos_folds: list of sets of pairs of positive interactions.
        :param pp_not_neg: set of low-quality positive interactions (non-negatives).
        :returns: list of sets of triples of the form (p1, p2, are_bound).
        """
        assert isinstance(pp_not_neg, set)

        folds = []
        for pos_fold in pos_folds:
            # Compute the candidate negative pairs
            fold_ps = set(p1 for p1, p2 in pos_fold) | \
                      set(p2 for p1, p2 in pos_fold)
            candidate_neg_pps = list(set(it.product(fold_ps, fold_ps)) - pp_not_neg)

            # Randomly sample a number of negative pairs equal to the number of
            # positives pairs
            pi = np.random.permutation(len(pos_fold))
            neg_fold = set(candidate_neg_pps[pi[i]] for i in xrange(len(fold_ps)))

            # Assemble the fold
            fold = set()
            fold.update((p1, p2, 1) for p1, p2 in pos_fold)
            fold.update((p1, p2, 0) for p1, p2 in neg_fold)
            folds.append(fold)

        return folds

    @staticmethod
    def _term_to_predicate(term):
        assert len(term.namespace)
        level = term.level if term.level >= 0 else None
        return "{}_{}_{}".format(term.namespace, term.level, term.id.replace(":", "_"))

    def _to_sbr_examples(self, ps, p_to_fun, labeled_pps):
        examples = []
        for p1, p2, state in labeled_pps:
            examples.append("ISPAIR({},{},{}{})".format(p1, p2, p1, p2))
            examples.append("ISPAIR({},{},{}{})".format(p1, p2, p2, p1))
            examples.append("BOUNDP({}{})={}".format(p1, p2, state))
        for p in ps:
            for term in p_to_fun[p]:
                examples.append("{}({})=1".format(self._term_to_predicate(term), p))
        return examples

    def _write_sbr_data(self, ps, p_to_fun, dag, folds):
        """Writes out the SBR data.

        .. warning::

            The interactions in the folds **must** be symmetric; this is not
            checked here.

        :param ps: list of protein IDs.
        :param p_to_fun: map between protein IDs and the ``GOTerm``'s they are annotated with.
        :param dag: the associated ``GODag``.
        :param folds: list of sets of triples of the form (protein ID, protein ID, are_bound).
        """
        # Write out the datapoints: proteins and protein pairs
        pos_pps, neg_pps = [], []
        for fold in folds:
            pos_pps.extend((p1, p2) for p1, p2, state in fold if state == 1)
            neg_pps.extend((p1, p2) for p1, p2, state in fold if state == 0)

        datapoints = []
        datapoints.extend("{}:PROTEIN".format(p) for p in ps)
        datapoints.extend("{}{}:PPAIR".format(p1, p2) for p1, p2 in (pos_pps + neg_pps))

        with open(os.path.join(self.dst, "sbr-datapoints"), "wb") as fp:
            fp.write("\n".join(datapoints))

        # Write out the predicates, which include (i) the BOUNDP predicate,
        # (ii) the ISPAIR predicate, and (iii) a predicate for each GO term.
        dag_terms = sorted(dag._id_to_term.values(), key = lambda term: term.name)

        predicates = []
        predicates.append("DEF ISPAIR(PROTEIN,PROTEIN,PPAIR);GIVEN;C;F")
        predicates.append("DEF BOUNDP(PPAIR);LEARN;C")
        predicates.extend("DEF {};LEARN;C".format(self._term_to_predicate(term))
                          for term in dag_terms)

        with open(os.path.join(self.dst, "sbr-predicates"), "wb") as fp:
            fp.write("\n".join(predicates))

        # Write out the rules
        # WRITEME

        # Write out the folds (examples)
        all_folds = set()
        for fold in folds:
            all_folds |= fold

        for k, fold in enumerate(folds):
            rest = all_folds - fold

            fold_ps = set(p1 for p1, _, _ in fold) | \
                      set(p2 for _, p2, _ in fold)

            rest_ps = set(p1 for p1, _, _ in rest) | \
                      set(p2 for _, p2, _ in rest)

            with open(os.path.join(self.dst, "sbr-fold{}-examples-train".format(k)), "wb") as fp:
                fp.write("\n".join(self._to_sbr_examples(rest_ps, p_to_fun, rest)))

            with open(os.path.join(self.dst, "sbr-fold{}-examples-test".format(k)), "wb") as fp:
                fp.write("\n".join(self._to_sbr_examples(fold_ps, p_to_fun, fold)))

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

    def run(self, min_sequence_len = 50, cdhit_threshold = 0.75):
        """Run the yeast prediction experiment."""

        # Retrieve the list of protein IDs (whose order will define the
        # structure of the kernels), and their mappings to feature IDs,
        # sequences and GO term annotations.
        ps          = self._cached(self._get_sgd_ids,
                                   "sgd_ids")
        p_to_feat   = self._cached(self._get_sgd_id_to_feat,
                                   "sgd_id_to_feat")
        p_to_seq    = self._cached(self._get_sgd_id_to_seq,
                                   "sgd_id_to_seq")
        p_to_term_ids = self._cached(self._get_sgd_id_to_fun,
                                   "sgd_id_to_fun")
        print _cls(self), ": found {} proteins".format(len(ps))

        for p in ps:
            assert p in p_to_seq, "'{}' has no sequence".format(p)
            assert p in p_to_term_ids, "'{}' has no GO annotation".format(p)
            assert p in p_to_feat, "'{}' has no associated feature ID".format(p)



        # Load the GO OBO file
        dag = GODag(os.path.join(self.src, "GO", "go-basic.obo"))

        # Fill in the GO data structure with the protein annotations and
        # propagate them to the root
        propagated_p_to_term_ids = self._cached(dag.annotate,
                                               "sgd_id_to_fun_propagated",
                                                p_to_term_ids,
                                                propagate = True)

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

        # Dump the dataset statistics prior to any preprocessing
        print _cls(self), ": dumping dataset statistics"
        self._dump_dataset_statistics(ps, p_to_seq, p_to_term_ids, dag, "raw")



        # Filter out sequences shorter than min_sequence_len residues, then
        # cluster the filtered sequences with CD-HIT
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= min_sequence_len, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), min_sequence_len)

        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, clusters = CDHit().run(filtered_p_seq, threshold = cdhit_threshold)
        filtered_ps = [list(cluster)[0][0] for cluster in clusters]
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(clusters), cdhit_threshold)

        # Preprocess the GO DAG
        filtered_p_to_term_ids = self._cached(dag.preprocess,
                                              "filtered_p_to_term_ids",
                                              filtered_ps,
                                              aspects = self._go_aspects,
                                              min_annot = self._min_go_annot,
                                              max_depth = self._max_go_depth)

        # Dump the dataset statistics after the preprocessing
        print _cls(self), ": dumping dataset statistics"
        self._dump_dataset_statistics(filtered_ps, p_to_seq, filtered_p_to_term_ids,
                                      dag, "preproc")



        # Query the high-quality protein-protein interactions, restricting
        # them to proteins in the filtered set
        pp_pos_hq, pp_pos_lq = self._get_sgd_pins(filtered_ps)

        # Create the folds based on function and positive interactions
        print _cls(self), ": computing the folds"
        pos_folds = self._cached(self._get_folds, "pos-folds",
                                 pp_pos_hq, filtered_p_to_term_ids)

        # Now that we have the function-balanced folds, pour some negative
        # interactions in
        print _cls(self), ": adding negatives"
        folds = self._cached(self._add_negatives, "folds", pos_folds, pp_pos_lq)



        # Write down the SBR input files
        # print _cls(self), ": writing SBR input files"
        # self._write_sbr_data(filtered_ps, filtered_p_to_term_ids, dag, folds)

        # Compute the protein kernels
        print _cls(self), ": computing the protein kernels"
        self._compute_protein_kernels(filtered_ps, p_to_feat, p_to_seq, folds)
