# -*- coding: utf-8 -*-

import numpy as np
from os.path import join
from itertools import product
from collections import defaultdict
from textwrap import dedent
from ocelot.services import _cls, CDHit
from ocelot.go import GODag
from ocelot.kernels import *
from ocelot.utils import permute
from ocelot import Experiment, Stage
from .yeast import *



class SGDExperiment(Experiment):
    """New experiment based on SGD and iPfam.

    This experiment is structured exactly the same as the Yip et al. [Yip09]_
    experiment, but performed on a different dataset based on SGD, BioGRID
    and iPfam.

    Parameters
    ----------
    go_aspects : list, optional
        List of GO aspects to retain, defaults to None (all of them).
    max_go_depth : int
        All GO annotations deeper than this are ignored, defaults to None
    min_go_annot : int
        All GO terms with fewer annotations than this are ignored, defaults to
        None
    All other parameters are passed to the parent class.
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

        print _cls(self), ": allowed GO ECs are '{}'".format(ECODES_TO_KEEP)

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

    def _load_sgd_resources(self):
        """Retrieves all SGD-derived resources.

        Returns
        -------
        ps : list
            Sorted list of protein (validated ORF) SGD IDs.
        p_to_feat : dict
            Map from SGD ID to protein (aka SGD feature) name.
        p_to_seq: dict
            Map from SGD ID to protein sequenze.
        p_to_term_ids : dict
            Map from SGD ID to GO term IDs.
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



    def _print_go_stats(self, dag, message):
        p_to_ids = dag.get_p_to_term_ids()

        id_to_ps = defaultdict(set)
        for p, ids in p_to_ids.iteritems():
            for id_ in ids:
                id_to_ps[id_].add(p)

        num_ps = len(p_to_ids)
        num_annot_ps = len(set(p for p, ids in p_to_ids.iteritems()
                               if len(ids) > 0))

        num_terms = len(id_to_ps)
        num_annot_terms = len(set(id_ for id_, ps in id_to_ps.iteritems()
                                  if len(ps) > 0))

        print _cls(self), ": {}/{} annotated proteins" \
                            .format(num_annot_ps, num_ps)
        print _cls(self), ": {}/{} annotated terms" \
                            .format(num_annot_terms, num_terms)

    def _get_go_annotations(self, filtered_ps, p_to_term_ids):
        """Load the GO OBO file, fill it in, and prune it."""

        print _cls(self), ": processing GO annotations (aspects={} max_depth={} min_annot={})" \
                            .format(self._go_aspects, self._min_go_annot,
                                    self._max_go_depth)

        # Keep only the filtered proteins
        p_to_term_ids = {p: term_ids for p, term_ids in p_to_term_ids.iteritems()
                         if p in filtered_ps}

        # Load the GO OBO file
        dag = GODag(join(self.src, "GO", "go.obo"))

        # Fill in the GO data structure with the protein annotations and
        # propagate them to the root
        dag.annotate(p_to_term_ids, propagate=True)
        self._print_go_stats(dag, "propagation")

        # Prune all useless terms from the GO
        dag.preprocess(filtered_ps, aspects=self._go_aspects,
                       min_annot=self._min_go_annot,
                       max_depth=self._max_go_depth)
        self._print_go_stats(dag, "preprocessing")

        return dag, dag.get_p_to_term_ids()


    def _check_ps_pps(self, ps, pps):
        """Checks that the protein pairs are (i) symmetric, and (ii) entirely
        contained in the list of proteins."""
        assert all((q, p) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"


    def _get_positive_pins(self, ps):
        """Computes the high-quality positive interactions, high+low-quality
        positive interactions, and the negative interactions."""

        # Query the high-quality protein-protein interactions
        pp_pos_hq = self._get_sgd_pin(ps=ps, manual_only=True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_ps_pps(ps, pp_pos_hq)

        # Query all (high+low-quality) protein-protein interactions
        pp_pos_lq = self._get_sgd_pin(ps=ps, manual_only=False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality PPIs (density = {})" \
                .format(len(pp_pos_lq), density)
        self._check_ps_pps(ps, pp_pos_lq)

        # Query (literally) all protein-protein actions annotated in STRING
        pp_pos_string = self._get_string_pin(ps=ps)
        density = float(len(pp_pos_string)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} STRING-quality PPIs (density = {})" \
                .format(len(pp_pos_string), density)
        self._check_ps_pps(ps, pp_pos_string)

        return pp_pos_hq, pp_pos_lq | pp_pos_string



    def _get_negative_pin(self, ps, pp_pos_hq, pp_pos_lq):
        """Computes the network of negative protein interactions.

        The idea is to sample the negative pairs from the complement of the
        positive PIN minus a set of candidate interactions (of any kind,
        really).
        """
        print _cls(self), ": sampling {} negative interactions".format(len(pp_pos_hq))
        pp_neg = set()
        while len(pp_neg) != len(pp_pos_hq):
            i, j = self._rng.randint(len(ps), size=2)
            p, q = ps[i], ps[j]
            if (p, q) not in pp_pos_lq and \
               (p, q) not in pp_pos_hq and \
               (p, q) not in pp_neg:
                pp_neg.update([(p, q), (q, p)])
        print _cls(self), ": got {} negative interactions".format(len(pp_neg))

        self._check_ps_pps(ps, pp_neg)
        return pp_neg,


    def _check_folds(self, ps, pps, folds):
        """Checks that the folds make sense."""

        # # All proteins must belong to at least a fold
        # ps_in_folds = set()
        # for fold in folds:
        #     ps_in_folds.update(p for p, _, _ in fold)
        #     ps_in_folds.update(q for _, q, _ in fold)
        # assert len(ps_in_folds) == len(set(ps))
        # del ps_in_folds

        # All protein pairs must belong to at least a fold
        pps_in_folds = set()
        for fold in folds:
            pps_in_folds.update((p, q) for p, q, _ in fold)
        assert len(pps_in_folds) == len(set(pps))
        del pps_in_folds

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

    def _compute_folds(self, ps, pp_pos, pp_neg, p_to_term_ids, num_folds=10):
        """Generates the folds.

        Folds are composed of *pairs* of (known interacting or assumed
        non-interacting) proteins.

        Protein pairs are split evenly according to their function.

        No two pairs can occur in distinct folds.

        However, individual proteins may occur in different folds, and some
        proteins may not occur altogether.

        Parameters
        ----------
        pp_pos : set
            Positive interactions.
        pp_neg : set
            Negative interactions.
        p_to_term_ids : dict
            Map between proteins and GO terms IDs.
        num_folds : int, optional
            Number of folds, defaults to 10.

        Returns
        -------
        folds : list
            Each fold is a set of triples of the form (protein, protein, state).
        """

        # XXX take into consideration INTRA interactions
        # XXX take into consideration unbound proteins
        # XXX take into consideration unannotated proteins

        def distribute_pps(folds, pps, state, p_to_term_ids):

            # Map from pairs to term IDs
            pp_to_term_ids = {}
            for p, q in pps:
                p_term_ids = p_to_term_ids.get(p, set())
                q_term_ids = p_to_term_ids.get(q, set())
                pp_to_term_ids[(p, q)] = p_term_ids | q_term_ids

            # Map from term IDs to pairs
            term_id_to_pps = defaultdict(set)
            for pp, term_ids in pp_to_term_ids.iteritems():
                if not len(term_ids):
                    term_id_to_pps["unannotated"].add(pp)
                else:
                    for term_id in term_ids:
                        term_id_to_pps[term_id].add(pp)

            # Distribute the interactions among the folds
            while len(term_id_to_pps):

                # Shuffle the folds
                folds = permute(folds, self._rng)

                # Pick the term with the least unallocated pairs
                cur_term_id, symm_pps = min(term_id_to_pps.iteritems(),
                                            key=lambda id_pps: len(id_pps[-1]))

                print "distributing interactions with term {} (# interactions = {})".format(cur_term_id, len(symm_pps))

                # Desymmetrize the pps to add
                asymm_pps = set()
                for p, q in symm_pps:
                    if not (q, p) in asymm_pps:
                        asymm_pps.add((p, q))

                # Evenly distribute the associated pairs among all folds.
                for i, (p, q) in enumerate(sorted(asymm_pps)):
                    folds[i % len(folds)].update([(p, q, state), (q, p, state)])

                # Update term_to_pps by removing all interactions that we
                # just allocated
                new_term_id_to_pps = {}
                for term_id in term_id_to_pps:
                    # Skip the term we just processed (it's done, right?)
                    if term_id == cur_term_id:
                        continue
                    new_pps = set(pp for pp in term_id_to_pps[term_id]
                                  if not pp in symm_pps)
                    if not len(new_pps):
                        continue
                    new_term_id_to_pps[term_id] = new_pps
                term_id_to_pps = new_term_id_to_pps

        print _cls(self), ": generating {} folds...".format(num_folds)
        folds = [set() for _ in range(num_folds)]

        print _cls(self), ": adding {} positive interactions...".format(len(pp_pos))
        distribute_pps(folds, pp_pos, True, p_to_term_ids)
        print _cls(self), ": adding {} negative interactions...".format(len(pp_neg))
        distribute_pps(folds, pp_neg, False, p_to_term_ids)

        # XXX proteins that do not interact (or non-interact) with anybody
        # do not appear in the folds.

        print _cls(self), ": checking fold sanity..."
        self._check_folds(ps, pp_pos | pp_neg, folds)

        return folds,



    def _dump_stats(self, ps, dag, p_to_term_ids, pp_pos_hq, pp_pos_lq, pp_neg,
                    folds):
        print "DATASET STATISTICS"

        print "# proteins =", len(ps)
        print "# hi-quality interactions =", len(pp_pos_hq)
        print "# lo-quality interactions =", len(pp_pos_lq)
        print "# negative interactions =", len(pp_neg)

        term_ids = set()
        for ids in p_to_term_ids.values():
            term_ids.update(ids)
        print "# annotated GO terms =", len(term_ids)

        print "GO DAG STATISTICS"

        level_to_term_ids = defaultdict(set)
        level_to_annots = defaultdict(list)
        level_to_ps = defaultdict(set)
        for id_ in term_ids:
            term = dag._id_to_term[id_]
            assert len(term.proteins)
            level_to_term_ids[term.level].add(id_)
            level_to_annots[term.level].extend(list(term.proteins))
            level_to_ps[term.level].update(term.proteins)

        for level in level_to_annots:
            num_terms = len(level_to_term_ids[level])
            num_bins = len([id_ for id_ in level_to_term_ids[level] if id_.startswith("BN")])
            num_annots = len(level_to_annots[level])
            num_ps = len(level_to_ps[level])
            print "level {} has {} terms (of which {} bins) with {} annotations over {} proteins".format(
                level, num_terms, num_bins, num_annots, num_ps)

        print "FOLD STATISTICS"

        p_to_folds = defaultdict(set)
        for k, fold in enumerate(folds):
            fold_ps = set()
            fold_ps.update(p for [p, _, _] in fold)
            fold_ps.update(q for [_, q, _] in fold)

            for p in fold_ps:
                p_to_folds[p].add(k)

            num_pos = len({(p, q) for p, q, s in fold if s})
            num_neg = len({(p, q) for p, q, s in fold if not s})

            print "# interactions in fold {}: {} over {} proteins (#pos={} #neg={} intra={} inter={})".format(
                k, len(fold), len(fold_ps), len(pos), len(neg))

        print "# proteins with GO annotations =", len(set(p for p in ps if p in p_to_term_ids and len(p_to_term_ids[p])))
        print "# proteins in a fold =", \
                len({p for p, folds in p_to_folds.items() if len(folds)})
        print "# proteins in at least 2 folds =", \
                len({p for p, folds in p_to_folds.items() if len(folds) >= 2})

        # XXX GO annotation cooccurrence

        # XXX GO annotation balance in folds

        # XXX draw full PPI (taken from the folds)

        # XXX draw folds over the full PPI

        return None,


    def _compute_p_colocalization_kernel(self, ps):
        p_to_context = self._get_sgd_id_to_context()
        contexts = [p_to_context[p] for p in ps]
        return self._compute_kernel(ColocalizationKernel, contexts, gamma=1.0)

    def _compute_p_gene_expression_kernel(self, ps, p_to_feat):
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}
        return self._compute_kernel(SGDGeneExpressionKernel,
                                    feat_to_i, self.src,
                                    ["Gasch_2000_PMID_11102521",
                                     "Spellman_1998_PMID_9843569"],
                                    do_normalize=True)

    def _compute_p_complex_kernel(self, ps, p_to_feat):
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}
        return self._compute_kernel(YeastProteinComplexKernel,
                                    feat_to_i, self.src)

    def _compute_p_interpro_kernel(self, ps):
        return self._compute_kernel(InterProKernel, ps,
                                    join(self.src, "interpro"), mode="match")

    def _compute_p_interpro_count_kernel(self, ps):
        return self._compute_kernel(InterProKernel, ps,
                                    join(self.src, "interpro"), mode="count")

    def _compute_p_pssm_kernel(self, ps, p_to_seq):
        return self._compute_kernel(PSSMKernel, ps, p_to_seq,
                                    join(self.src, "pssm"), k=4, threshold=6.0)

    def _compute_average_kernel(self, *matrices):
        return sum(matrices) / len(matrices),



    def _compute_p_pp_order(self, ps, pp_pos, pp_neg):
        p_to_i = {p: i for i, p in enumerate(sorted(ps))}

        pps = set()
        pps.update(pp_pos)
        pps.update(pp_neg)
        pp_indices = sorted([(p_to_i[p], p_to_i[q]) for p, q in sorted(pps)])

        return p_to_i, pp_indices

    def _compute_pp_kernel(self, pp_indices, submatrix):
        return self._compute_kernel(PairwiseKernel, pp_indices, submatrix),


    def _write_sbr_datapoints(self, p_to_i, pp_indices):
        i_to_p = {i: p for p, i in p_to_i.iteritems()}

        lines = []
        for p, i in sorted(p_to_i.items(), key=lambda p_i: p_i[-1]):
            lines.append("{};PROTEIN".format(p))
        for i, j in pp_indices:
            lines.append("{}-{};PPAIR".format(i_to_p[i], i_to_p[j]))

        with open(join(self.dst, "sbr-datapoints.txt"), "wb") as fp:
            fp.write("\n".join(lines))

    def _write_sbr_predicates(self, dag):
        lines = []
        lines.append("DEF BOUND(PPAIR);LEARN;C")
        lines.append("DEF IN(PROTEIN,PROTEIN,PPAIR);GIVEN;C;F")
        lines.extend("DEF {};LEARN;C".format(term_to_predicate(term))
                     for term in dag._id_to_term.itervalues())

        with open(join(self.dst, "sbr-predicates.txt"), "wb") as fp:
            fp.write("\n".join(lines))

    def _write_sbr_rules(self, dag):
        for aspect in dag.ASPECTS:

            # Parents imply the OR of the children
            lines = []
            for term in dag._id_to_term.itervalues():
                if term.namespace != aspect:
                    continue
                children = [child for child, _ in term.get_children(dag)]
                if not len(children):
                    continue
                children_or = " OR ".join(term_to_predicate(child) + "(p)"
                                          for child in children)
                lines.append("forall p [({}(p)) => ({})];LEARN;L1;LUKASIEWICZ_TNORM;1;1" \
                                .format(term_to_predicate(term), children_or))

            with open(join(self.dst, "sbr-rules-{}-term_implies_or_of_children.txt".format(aspect)), "wb") as fp:
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
                lines.append("forall p [({}(p)) => ({})];LEARN;L1;LUKASIEWICZ_TNORM;1;1" \
                                .format(term_to_predicate(term), parents_and))
            with open(join(self.dst, "sbr-rules-{}-term_implies_and_of_parents.txt".format(aspect)), "wb") as fp:
                fp.write("\n".join(lines))

            # Interaction implies same function within the same level
            aspect_terms_by_level = defaultdict(set)
            for term in dag.get_terms_by_level():
                if term.namespace == aspect:
                    aspect_terms_by_level[term.level].add(term)

            lines = []
            for level in sorted(aspect_terms_by_level):
                head = " OR ".join("({}(p) AND {}(q))" \
                                    .format(term_to_predicate(term),
                                            term_to_predicate(term))
                                   for term in aspect_terms_by_level[level])
                lines.append("forall p forall q forall pq [(BOUND(pq)) AND (IN(p,q,pq)) => {}];LEARN;L1;MINIMUM_TNORM;1;1;IN".format(head))
            with open(join(self.dst, "sbr-rules-interaction_implies_same_{}_levelwise.txt".format(aspect)), "wb") as fp:
                fp.write("\n".join(lines))

    def _write_sbr_dataset(self, ps, dag, p_to_term_ids, pp_pos, pp_neg,
                           p_to_i, pp_indices, folds):
        """Writes the SBR dataset."""

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

        self._write_sbr_datapoints(p_to_i, pp_indices)
        self._write_sbr_predicates(dag)
        self._write_sbr_rules(dag)

        def fold_to_ps(fold):
            return set(p for p, _, _ in fold) | set(q for _, q, _ in fold)

        def write_functions(path, ps, dag, p_to_term_ids):
            """Writes GO term predicates for each protein in the fold."""
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
            """Writes IN/3 and BOUND/2 predicates for all (p, q, state) tuples
            in the fold."""
            lines = []
            for p, q, state in fold:
                lines.append("IN({},{},{}-{})=1".format(p, q, p, q))
                lines.append("BOUND({}-{})={}".format(p, q, {True:1, False:0}[state]))
            with open(path, "wb") as fp:
                fp.write("\n".join(lines))

        p_to_term_ids = dag.get_p_to_term_ids()

        # Write the examples for each fold
        for k, test_set in enumerate(folds):
            l = (k + 1) % len(folds)
            validation_set = folds[l]

            train_set = set()
            for i, fold in enumerate(folds):
                if i != k and i != l:
                    train_set.update(fold)

            ps_in_test_set = fold_to_ps(test_set)
            ps_in_validation_set = fold_to_ps(validation_set) - ps_in_test_set
            ps_in_train_set = fold_to_ps(train_set) - \
                (ps_in_test_set | ps_in_validation_set)

            print _cls(self), \
                "fold {} examples: #ps = {}, {}, {}; #pps = {}, {}, {}" \
                    .format(k, len(ps_in_train_set),
                            len(ps_in_validation_set), len(ps_in_test_set),
                            len(train_set), len(validation_set),
                            len(test_set))

            write_functions(join(self.dst, "sbr-fold{}-testset-fun.txt".format(k)),
                            ps_in_test_set, dag, p_to_term_ids)
            write_interactions(join(self.dst, "sbr-fold{}-testset-int.txt".format(k)),
                               test_set)

            write_functions(join(self.dst, "sbr-fold{}-validset-fun.txt".format(k)),
                            ps_in_validation_set, dag, p_to_term_ids)
            write_interactions(join(self.dst, "sbr-fold{}-validset-int.txt".format(k)),
                               validation_set)

            write_functions(join(self.dst, "sbr-fold{}-trainset-fun.txt".format(k)),
                            ps_in_train_set, dag, p_to_term_ids)
            write_interactions(join(self.dst, "sbr-fold{}-trainset-int.txt".format(k)),
                               train_set)

        return True,



    def run(self, min_sequence_len=50, cdhit_threshold=0.75):
        """Run the yeast prediction experiment.

        Parameters
        ----------
        min_sequence_len : int, optional
            Minimum allowed sequence length, defaults to 50.
        cdhit_threshold : float, optional
            CD-HIT sequence similarity threshold, defaults to 0.75.
        """
        STAGES = (

            Stage(self._load_sgd_resources,
                  [], ['ps', 'p_to_feat', 'p_to_seq', 'p_to_term_ids']),

            Stage(self._filter_ps,
                  ['ps', 'p_to_seq', 'min_sequence_len', 'cdhit_threshold'],
                  ['filtered_ps']),

            Stage(self._get_go_annotations,
                  ['filtered_ps', 'p_to_term_ids'],
                  ['filtered_dag', 'filtered_p_to_term_ids']),

            Stage(self._get_positive_pins,
                  ['filtered_ps'], ['pp_pos_hq', 'pp_pos_lq']),

            Stage(self._get_negative_pin,
                  ['filtered_ps', 'pp_pos_hq', 'pp_pos_lq'],
                  ['pp_neg']),

            Stage(self._compute_folds,
                  ['filtered_ps', 'pp_pos_hq', 'pp_neg', 'filtered_p_to_term_ids'],
                  ['folds']),

            Stage(self._dump_stats,
                  ['filtered_ps', 'filtered_dag', 'filtered_p_to_term_ids',
                   'pp_pos_hq', 'pp_pos_lq', 'pp_neg', 'folds'],
                  ['__dummy_stats']),

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

            Stage(self._compute_average_kernel,
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
                  ['__dummy_p_kernels']),

            Stage(self._compute_p_pp_order,
                  ['filtered_ps', 'pp_pos_hq', 'pp_neg'],
                  ['filtered_p_to_i', 'pp_indices']),

            Stage(self._compute_pp_kernel,
                  ['filtered_ps', 'pp_indices', 'p_average_kernel'],
                  ['pp_average_kernel']),

            Stage(lambda *args, **kwargs: None,
                  [
                    'pp_average_kernel',
                  ],
                  ['__dummy_pp_kernels']),

            Stage(self._write_sbr_dataset,
                  ['filtered_ps', 'filtered_dag', 'filtered_p_to_term_ids',
                   'pp_pos_hq', 'pp_neg', 'pp_indices', 'folds'],
                  ['__dummy_sbr']),
        )

        TARGETS = ('__dummy_stats', '__dummy_p_kernels', '__dummy_sbr')

        context = {
            "min_sequence_len": min_sequence_len,
            "cdhit_threshold": cdhit_threshold,
        }
        self._scheduler.run(STAGES, TARGETS, context=context)
