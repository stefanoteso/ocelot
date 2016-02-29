# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import product
from collections import defaultdict
from ocelot.services import _cls, CDHit
from ocelot.go import GODag, GOTerm
from ocelot.kernels import *
from ocelot.scheduler import Stage
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
                  for bindings in self.iterquery(query, n=1))
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
        for bindings in self.iterquery(query, n=2):
            sgd_id = bindings[u"orf"].split(".")[-1]
            sgd_id_to_seq[sgd_id] = bindings[u"seq"]
        return sgd_id_to_seq

    def _get_sgd_id_to_term_ids(self):
        ECODES_TO_KEEP = [
            # Experimental evidence codes. Use of an experimental evidence code
            # in a GO annotation indicates that the cited paper displayed
            # results from a physical characterization of a gene or gene
            # product that has supported the association of a GO term.
            "EXP",  # Inferred from Experiment (EXP)
            "IDA",  # Inferred from Direct Assay (IDA)
            "IPI",  # Inferred from Physical Interaction (IPI)
            "IMP",  # Inferred from Mutant Phenotype (IMP)
            "IGI",  # Inferred from Genetic Interaction (IGI)
            "IEP",  # Inferred from Expression Pattern (IEP)

            # Computational analysis evidence codes. Use of the computational
            # analysis evidence codes indicates that the annotation is based on
            # an in silico analysis of the gene sequence and/or other data as
            # described in the cited reference.  The evidence codes in this
            # category also indicate a varying degree of curatorial input.
        #    "ISS",  # Inferred from Sequence or structural Similarity (ISS)
        #    "ISO",  # Inferred from Sequence Orthology (ISO)
        #    "ISA",  # Inferred from Sequence Alignment (ISA)
        #    "ISM",  # Inferred from Sequence Model (ISM)
        #    "IGC",  # Inferred from Genomic Context (IGC)
        #    "IBA",  # Inferred from Biological aspect of Ancestor (IBA)
        #    "IBD",  # Inferred from Biological aspect of Descendant (IBD)
        #    "IKR",  # Inferred from Key Residues (IKR)
        #    "IRD",  # Inferred from Rapid Divergence(IRD)
        #    "RCA",  # Inferred from Reviewed Computational Analysis (RCA)

            # Author statement codes. Author statement codes indicate that the
            # annotation was made on the basis of a statement made by the
            # author(s) in the reference cited.
            "TAS",  # Traceable Author Statement (TAS)
        #    "NAS",  # Non-traceable Author Statement (NAS)

            # Curatorial evidence codes. Use of the curatorial statement
            # evidence codes indicates an annotation made on the basis of a
            # curatorial judgement that does not fit into one of the other
            # evidence code classifications.
            "IC",   # Inferred by Curator (IC)
        #    "ND",   # No biological Data available (ND) evidence code

            # Automatically-assigned evidence code. Assigned by automated
            # methods, without curatorial judgement
        #    "IEA",  # Inferred from Electronic Annotation (IEA)
        ]

        query = """
        SELECT ?orf ?fun ?ecode
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_qualifier ocelot:sgd_id_qualifier.Verified .
            ?annot a ocelot:sgd_go_annotation ;
                ocelot:sgd_go_annotation_has_ecode ?ecode ;
                ocelot:sgd_go_annotation_has_sgd_id ?orf ;
                ocelot:sgd_go_annotation_has_term ?fun .
        }}
        """
        sgd_id_to_fun = defaultdict(set)
        for bindings in self.iterquery(query, n=3):
            sgd_id  = bindings[u"orf"].split(".")[-1]
            fun     = bindings[u"fun"].split("#")[-1]
            ecode   = bindings[u"ecode"]
            if str(ecode) in ECODES_TO_KEEP:
                sgd_id_to_fun[sgd_id].add(fun)
        return dict(sgd_id_to_fun)

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
        for bindings in self.iterquery(query, n=2):
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
        for bindings in self.iterquery(query, n=5):
            orf     = bindings[u"orf"].split(".")[-1]
            chrom   = bindings[u"chrom"].split(".")[-1]
            strand  = bindings[u"strand"]
            start   = bindings[u"start"]
            stop    = bindings[u"stop"]
            assert strand in ("C", "W")
            p_to_context[orf] = \
                (chrom, min(start, stop) + 0.5 * np.fabs(start - stop))
        return p_to_context

    def _get_sgd_pin(self, ps=None, manual_only=False):
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
        for bindings in self.iterquery(query, n=2):
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

    def _get_string_pin(self, ps=None):
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
        for bindings in self.iterquery(query, n=2):
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



    def _get_sgd_pins(self, ps):
        """Computes the high-quality positive interactions, high+low-quality
        positive interactions, and the negative interactions."""

        def check_ps_pps(ps, pps):
            """Checks that the protein pairs are (i) symmetric, and (ii) entirely
            contained in the list of proteins."""
            assert all((q, p) in pps for (p, q) in pps), \
                "pairs are not symmetric"
            assert all(p in ps for p, _ in pps), \
                "singletons and pairs do not match"

        pp_pos_hq = self._get_sgd_pin(ps=ps, manual_only=True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        check_ps_pps(ps, pp_pos_hq)

        # Query all (high+low-quality) protein-protein interactions
        pp_pos_lq = self._get_sgd_pin(ps=ps, manual_only=False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality PPIs (density = {})" \
                .format(len(pp_pos_lq), density)
        check_ps_pps(ps, pp_pos_lq)

        # Query (literally) all protein-protein actions annotated in STRING
        pp_pos_string = self._get_string_pin(ps=ps)
        density = float(len(pp_pos_string)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} STRING-quality PPIs (density = {})" \
                .format(len(pp_pos_string), density)
        check_ps_pps(ps, pp_pos_string)

        return pp_pos_hq, pp_pos_lq | pp_pos_string



    def _compute_negative_pin(self, ps, pp_pos_hq, pp_pos_lq):
        """Computes the network of negative protein interactions.

        The idea is to sample the negative pairs from the complement of the
        positive PIN minus a set of candidate interactions (of any kind,
        really).
        """

        def symmetrize(pairs):
            return set(pairs) | set((q, p) for p, q in pairs)

        def desymmetrize(pairs):
            pairs_asymm = set()
            for p, q in pairs_asymm:
                if not (q, p) in pairs:
                    pairs_asymm.add((p, q))
            return pairs_asymm

        # Take the complement of the symmetrical positives, subtract the known
        # positives (both high-quality and low-quality, just to be sure); the
        # complement will be symmetrical, since the HQ and LQ interactions are
        print _cls(self), ": computing the complement of the positive PIN"
        pp_neg_all = set(product(ps, ps)) - (pp_pos_lq | pp_pos_hq)

        print _cls(self), ": computing asymmetric negative PIN"
        pp_neg_asymm = list(desymmetrize(pp_neg_all))

        print _cls(self), ": computing asymmetric positive PIN"
        pp_pos_asymm = desymmetrize(pp_pos_hq)

        del pp_pos_hq
        del pp_pos_lq
        del pp_neg_all

        # Sample the negative interactions from the complement
        print _cls(self), ": sampling the negative PIN"
        pi = self._rng.permutation(len(pp_neg_asymm))
        pp_neg_half = set(pp_neg_asymm[pi[i]] for i in xrange(len(pp_pos_asymm)))

        # Symmetrize the negatives
        pp_neg = symmetrize(pp_neg_half)
        self._check_p_pp_are_sane(ps, pp_neg)

        return pp_neg,



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
                                    mode="match")

    def _compute_p_interpro_count_kernel(self, ps):
        return self._compute_kernel(InterProHistogramKernel, ps, self.dst,
                                    mode="count")

    def _compute_p_pssm_kernel(self, ps, p_to_seq):
        return self._compute_kernel(PSSMKernel, ps, p_to_seq, self.dst,
                                    k=4, threshold=6.0)

    def _compute_p_average_kernel(self, *matrices):
        return sum(matrices) / len(matrices)

    def _compute_pp_kernel(self, ps, folds, submatrix):
        p_to_i = {p: i for i, p in enumerate(ps)}

        pps = set()
        for fold in folds:
            pps.update((p1, p2) for p1, p2, _ in fold)
        pps = sorted(pps)

        pp_indices = [(p_to_i[p1], p_to_i[p2]) for p1, p2 in pps]
        return self._compute_kernel(PairwiseKernel, pp_indices, submatrix)



    def _compute_folds(self, pp_pos_hq, pp_pos_lq, p_to_term_ids, num_folds=10):
        """Generates the folds.

        The interactions are split into folds according to their functions.
        No two protein-protein interactions can occur in two distinct folds.

        Parameters
        ----------
        pp_pos_hq : set
            High-quality positive interactions.
        pp_pos_lq : set
            Low-quality positive interactions.
        p_to_term_ids : dict
            Map between proteins and GO terms IDs.
        num_folds : int, optional
            Number of folds, defaults to 10.

        Returns
        -------
        folds : list
            Each fold is a set of triples of the form (protein, protein, state).
        """

        def permute(l, rng):
            pi = list(rng.permutation(len(l))
            return [l[pi[i]] for i in xrange(len(l))]

        def check_folds(folds):

            # Interactions must be symmetric
            for k, fold in enumerate(folds):
                assert all((q, p, state) in fold for (p, q, state) in fold), \
                    "fold {} is not symmetric".format(k)

            # Folds must not overlap
            for (k1, fold1), (k2, fold2) in product(enumerate(folds), enumerate(folds)):
                if k1 >= k2:
                    continue
                for (p1, q1, state1), (p2, q2, state) in product(fold1, fold2):
                    assert (p1, q1) != (p2, q2), \
                        "folds {} and {} overlap".format(k1, k2)

            # Interactions must be either positive or negative, but not both
            interactions = set()
            for fold in folds:
                interactions.update(fold)
            for p, q, state in interactions:
                assert not(p, q, not state) in interactions, \
                    "folds are inconsistent"

        # Map from high-quality interacting protein pairs to the GO terms that
        # either protein is annotated with, and vice-versa
        pp_to_term_ids = {(p1, p2): set(p_to_term_ids[p1]) | set(p_to_term_ids[p2])
                       for p1, p2 in pp_pos_hq}

        # Map from each GO term ID to the protein pairs it is annotated with
        term_id_to_pps = defaultdict(set)
        for pp, term_ids in pp_to_term_ids.iteritems():
            for term_id in term_ids:
                term_id_to_pps[term_id].add(pp)

        # Assign an index to all term IDs
        term_id_to_i = {term_id: i for i, term_id in enumerate(term_id_to_pps.keys())}

        # Compute the ideal (perfect average) GO term distribution
        average = np.zeros((len(term_id_to_pps),))
        for pp, term_ids in pp_to_term_ids.iteritems():
            for term_id in term_ids:
                average[term_id_to_i[term_id]] += 1
        average *= 1.0 / num_folds

        # Sanity check. Note that self-interactions are fine.
        for term_id, term_pps in term_id_to_pps.iteritems():
            for p1, p2 in term_pps:
                assert (p2, p1) in term_pps

        # Generate the folds
        print _cls(self), ": generating the folds..."

        folds = [set() for _ in range(num_folds)]
        while len(term_id_to_pps):

            # Pick the term with the least unprocessed protein-protein pairs
            cur_term_id, cur_pps = min(term_id_to_pps.iteritems(),
                                       key=lambda term_id_and_pps: len(term_id_and_pps[1]))
            print _cls(self), ": best term: {} (num pairs = {})".format(cur_term_id, len(cur_pps))

            # Evenly distribute the associated protein-protein pairs among all
            # folds, taking into account the fact that if (p1, p2) is in
            # cur_pps, then (p2, p1) is in cur_pps as well. Also, make sure
            # to allow self-interactions.
            pps_to_add = set()
            for p1, p2 in cur_pps:
                if not (p2, p1) in pps_to_add:
                    pps_to_add.add((p1, p2))
            folds = permute(folds, self._rng)
            for i, (p1, p2) in enumerate(pps_to_add):
                folds[i % num_folds].update([(p1, p2, True), (p2, p1, True)])

            # Update term_to_pps
            new_term_id_to_pps = {}
            for term_id in term_id_to_pps:
                # Skip the term we just processed
                if term_id == cur_term_id:
                    continue
                # Skip the protein pairs annotated with it
                new_term_pps = set(pp for pp in term_id_to_pps[term_id]
                                   if not pp in cur_pps)
                if len(new_term_pps) == 0:
                    continue
                new_term_id_to_pps[term_id] = new_term_pps
            term_id_to_pps = new_term_id_to_pps

        print _cls(self), ": checking fold sanity..."
        check_folds(folds)

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
            candidate_neg_pps = set(product(ps_in_fold, ps_in_fold)) - candidate_pos_pps - all_candidate_neg_pps

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
        check_folds(folds)

        return folds,

    def _write_sbr_dataset(self, ps, dag, p_to_term_ids, folds):
        """Writes the SBR dataset."""

        # Figure out how the pairs are laid out in the pairwise kernels
        p_to_i = {p: i for i, p in enumerate(ps)}

        pps = set()
        for fold in folds:
            pps.update((p1, p2) for p1, p2, _ in fold)
        pps = sorted(pps)

        def term_to_predicate(term):
            assert len(term.namespace)
            assert term.level >= 0
            d = {
                "namespace": term.namespace,
                "level": term.level,
                "id": term.id.replace(":", "-"),
                "name": term.name.replace(" ", "-").replace(":", "-").replace(",", "")
            }
            return "{namespace}-lvl{level}-{id}-{name}".format(**d)

        # Write the datapoints, including proteins and protein pairs
        lines = []
        lines.extend("{};PROTEIN;1:1".format(p) for p in ps)
        lines.extend("{}-{};PPAIR;1:1".format(pp[0], pp[1]) for pp in pps)
        with open(os.path.join(self.dst, "sbr-datapoints.txt"), "wb") as fp:
            fp.write("\n".join(lines))

        # Write the predicates, including BOUNDP, ISPAIR and per-term predicates
        lines = []
        lines.extend("DEF {};LEARN;C".format(term_to_predicate(term))
                     for term_id, term in dag._id_to_term.iteritems())
        lines.append("DEF BOUND(PPAIR);LEARN;C")
        lines.append("DEF IN(PROTEIN,PROTEIN,PPAIR);GIVEN;C;F")
        with open(os.path.join(self.dst, "sbr-predicates.txt"), "wb") as fp:
            fp.write("\n".join(lines))

        # Write the rules
        for aspect in ("biological_process", "cellular_component", "molecular_function"):

            # Parents imply the OR of the children
            lines = []
            for term in dag._id_to_term.itervalues():
                if term.namespace != aspect:
                    continue
                children = [child for child, _ in term.get_children(dag)]
                if not len(children):
                    continue
                children_or = " OR ".join(term_to_predicate(child) + "(p)" for child in children)
                lines.append("forall p [({}(p)) => ({})];LEARN;L1;LUKASIEWICZ_TNORM;1;1".format(term_to_predicate(term), children_or))
            with open(os.path.join(self.dst, "sbr-rules-{}-term_implies_or_of_children.txt".format(aspect)), "wb") as fp:
                fp.write("\n".join(lines))

            # Children imply their parents
            lines = []
            for term in dag._id_to_term.itervalues():
                if term.namespace != aspect:
                    continue
                parents = [parent for parent, _ in term.get_parents(dag)]
                if not len(parents):
                    continue
                parents_and = " AND ".join(term_to_predicate(parent) + "(p)" for parent in parents)
                lines.append("forall p [({}(p)) => ({})];LEARN;L1;LUKASIEWICZ_TNORM;1;1".format(term_to_predicate(term), parents_and))
            with open(os.path.join(self.dst, "sbr-rules-{}-term_implies_and_of_parents.txt".format(aspect)), "wb") as fp:
                fp.write("\n".join(lines))

            # Interaction implies same function within the same level
            aspect_terms_by_level = defaultdict(set)
            for term in dag.get_terms_by_level():
                if term.namespace == aspect:
                    aspect_terms_by_level[term.level].add(term)

            lines = []
            for level in sorted(aspect_terms_by_level):
                head = " OR ".join("({}(p) AND {}(q))".format(term_to_predicate(term), term_to_predicate(term)) for term in aspect_terms_by_level[level])
                lines.append("forall p forall q forall pq [(BOUND(pq)) AND (IN(p,q,pq)) => {}];LEARN;L1;MINIMUM_TNORM;1;1;IN".format(head))
            with open(os.path.join(self.dst, "sbr-rules-interaction_implies_same_{}_levelwise.txt".format(aspect)), "wb") as fp:
                fp.write("\n".join(lines))

        def fold_to_ps(fold):
            return set(p for p, _, _ in fold) | set(q for _, q, _ in fold)

        def write_functions(path, ps, dag, p_to_term_ids):
            lines = []
            for p in ps:
                if not p in p_to_term_ids:
                    continue
                assert len(p_to_term_ids) > 0
                for term in dag._id_to_term.itervalues():
                    state = {True: "1", False: "0"}[term.id in p_to_term_ids[p]]
                    lines.append("{}({})={}".format(term_to_predicate(term), p, state))
            with open(path, "wb") as fp:
                fp.write("\n".join(lines))

        def write_interactions(path, fold):
            lines = []
            for p, q, state in fold:
                lines.append("IN({},{},{}-{})=1".format(p, q, p, q))
                lines.append("BOUND({}-{})={}".format(p, q, state))
            with open(path, "wb") as fp:
                fp.write("\n".join(lines))

        p_to_term_ids = dag.get_p_to_term_ids()

        # Write the examples for each fold
        for k, test_fold in enumerate(folds):
            l = (k + 1) % len(folds)

            validation_fold = folds[l]

            training_set = set()
            for i, fold in enumerate(folds):
                if i != k and i != l:
                    training_set.update(fold)

            ps_in_training_set = fold_to_ps(training_set)
            ps_in_test_set = fold_to_ps(test_fold) - ps_in_training_set
            ps_in_validation_set = fold_to_ps(validation_fold) - (ps_in_training_set | ps_in_test_set)

            print _cls(self), \
                "fold {} examples: #ps = {}, {}, {}; #pps = {}, {}, {}" \
                    .format(k, len(ps_in_training_set),
                            len(ps_in_validation_set), len(ps_in_test_set),
                            len(training_set), len(validation_fold),
                            len(test_fold))

            write_functions(os.path.join(self.dst, "sbr-fold{}-testset-fun.txt".format(k)),
                            ps_in_test_set, dag, p_to_term_ids)
            write_functions(os.path.join(self.dst, "sbr-fold{}-validset-fun.txt".format(k)),
                            ps_in_validation_set, dag, p_to_term_ids)
            write_functions(os.path.join(self.dst, "sbr-fold{}-trainset-fun.txt".format(k)),
                            ps_in_training_set, dag, p_to_term_ids)

            write_interactions(os.path.join(self.dst, "sbr-fold{}-testset-int.txt".format(k)),
                               test_fold)
            write_interactions(os.path.join(self.dst, "sbr-fold{}-validset-int.txt".format(k)),
                               validation_fold)
            write_interactions(os.path.join(self.dst, "sbr-fold{}-trainset-int.txt".format(k)),
                               training_set)

        return True,

    def _load_sgd_resources(self):
        """Retrieves all SGD-derived resources.

        Namely:

        - the list of protein (validated ORF) IDs.
        - the ID -> yeast protein name map.
        - the ID -> sequence map.
        - the ID -> GO term annotations map.

        :returns: the four items above.
        """
        def to_link(p):
            return "http://www.yeastgenome.org/locus/{}/go".format(p)

        ps = self._get_sgd_ids()
        print _cls(self), ": found {} proteins (validated ORFs)".format(len(ps))

        p_to_feat = self._get_sgd_id_to_feat()
        print _cls(self), ": found {} proteins with known SGD feature".format(len(p_to_feat))

        p_to_seq = self._get_sgd_id_to_seq()
        print _cls(self), ": found {} proteins with known SGD sequence".format(len(p_to_seq))

        p_to_term_ids = self._get_sgd_id_to_term_ids()
        print _cls(self), ": found {} proteins with known GO terms".format(len(p_to_term_ids))

        for p in ps:
            assert p in p_to_seq, "'{}' has no sequence".format(p)
            assert p in p_to_feat, "'{}' has no associated feature ID".format(p)
            if not p in p_to_term_ids: print "'{}' has no associated GO terms".format(to_link(p))

        return ps, p_to_feat, p_to_seq, p_to_term_ids

    def _filter_ps(self, ps, p_to_seq, min_sequence_len, cdhit_threshold):
        """Filter out short sequences and cluster them with CD-HIT."""

        # Filter out sequences shorter than min_sequence_len residues
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= min_sequence_len, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), min_sequence_len)

        # Cluster the filtered sequences with CD-HIT
        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, clusters = CDHit().run(filtered_p_seq, threshold=cdhit_threshold)
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(clusters), cdhit_threshold)

        # Take one protein from each cluster
        filtered_ps = [list(cluster)[0][0] for cluster in clusters]
        print _cls(self), ": there are {} filtered proteins".format(len(filtered_ps))

        return sorted(filtered_ps),

    def _load_go_and_filter(self, filtered_ps, p_to_term_ids):
        """Load the GO OBO file, fill it in, and prune it."""

        # Load the GO OBO file
        dag = GODag(os.path.join(self.src, "GO", "go.obo"))

        # Fill in the GO data structure with the protein annotations and
        # propagate them to the root
        dag.annotate(p_to_term_ids, propagate=True)

        id_to_term = dag.get_id_to_term()
        print _cls(self), ": '{}' GO annotations after propagation" \
                            .format(sum(len(id_to_term[id_].proteins)
                                        for id_ in id_to_term.iterkeys()))

        # Prune all useless terms from the GO
        print _cls(self), ": preprocessing the DAG: aspects={} max_depth={} min_annot={}" \
                            .format(self._go_aspects, self._min_go_annot,
                                    self._max_go_depth)
        dag.preprocess(filtered_ps, aspects=self._go_aspects,
                       min_annot=self._min_go_annot,
                       max_depth=self._max_go_depth)

        id_to_term = dag.get_id_to_term()
        print _cls(self), ": '{}' GO annotations after preprocessing" \
                            .format(sum(len(id_to_term[id_].proteins)
                                        for id_ in id_to_term.iterkeys()))

        # At this point, some proteins may not be annotated anymore. Fill in
        # the blanks by mapping them to an empty set.
        filtered_p_to_term_ids = dict(dag.get_p_to_term_ids())
        for p in set(filtered_ps) - set(filtered_p_to_term_ids.keys()):
            filtered_p_to_term_ids[p] = set()

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

            Stage(self._compute_p_interpro_count_kernel,
                  ['filtered_ps'], ['p_interpro_count_kernel']),

            Stage(self._compute_p_pssm_kernel,
                  ['filtered_ps', 'p_to_seq'], ['p_pssm_kernel']),

            Stage(self._compute_p_average_kernel,
                  [
                    'p_colocalization_kernel',
                    'p_gene_expression_kernel',
                    'p_complex_kernel',
                    'p_interpro_kernel',
                    'p_pssm_kernel'
                  ],
                  ['p_average_kernel']),

            Stage(lambda *args, **kwargs: None,
                  [
                    'p_colocalization_kernel',
                    'p_gene_expression_kernel',
                    'p_complex_kernel',
                    'p_interpro_kernel',
                    'p_interpro_count_kernel',
                    'p_pssm_kernel',
                    'p_average_kernel',
                  ],
                  ['p_kernels']),

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

            Stage(self._write_sbr_dataset,
                  ['filtered_ps', 'filtered_dag', 'filtered_p_to_term_ids', 'folds'],
                  ['sbr_dataset']),
        )

        TARGETS = (
            'sbr_dataset',
        )

        context = {
            "min_sequence_len": min_sequence_len,
            "cdhit_threshold": cdhit_threshold,
            "dump_states": dump_stats,
        }
        self._scheduler.run(STAGES, TARGETS, context=context)
