# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
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
        super(SGDExperiment, self).__init__(*args, **kwargs)

    def _get_sgd_ids(self):
        query = """
        SELECT ?orf ?seq
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF .
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
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF .
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
        assert all((p, q) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"

    def _get_sgd_pins(self, ps):
        """Computes the high-quality positive interactions, high+low-quality
        positive interactions, and the negative interactions."""
        pp_pos_hq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_hq.pickle",
                                 ps = ps, manual_only = True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_p_pp_are_sane(ps, pp_pos_hq)

        # Query all (high+low-quality) protein-protein interactions
        pp_pos_lq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_lq.pickle",
                                 ps = ps, manual_only = False)
        density = float(len(pp_pos_lq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} lo-quality PPIs (density = {})" \
                .format(len(pp_pos_lq), density)
        self._check_p_pp_are_sane(ps, pp_pos_lq)

        # Here we pass the low-quality interactions so as to get a better
        # approximation of the negative set.
        pp_neg = self._cached(self._get_negative_pin,
                              "sgd_id_interactions_neg.pickle",
                              ps, pp_pos_lq)
        print _cls(self), ": sampled {} p-p negative interactions" \
                .format(len(pp_neg))
        self._check_p_pp_are_sane(ps, pp_neg)

        return pp_pos_hq, pp_pos_lq, pp_neg

    def _get_sgd_din(self):
        query = """
        SELECT DISTINCT ?id ?family ?chain ?evalue ?complex
        FROM <{default_graph}>
        WHERE {{
            # Pick all SGD ORFs
            ?orf a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
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

    def _compute_protein_kernels(self, ps, p_to_feat):
        feat_to_i = {p_to_feat[p]: i for i, p in enumerate(ps)}

        # Compute the gene colocalization kernel
        p_to_context = self._cached(self._get_sgd_id_to_context,
                                    "sgd_id_to_context.pickle")
        contexts = [p_to_context[p] for p in ps]
        self._cached_kernel(ColocalizationKernel, len(ps),
                            "p-colocalization-kernel",
                            contexts, gamma = 1.0)

        # Compute the gene expression kernel
        self._cached_kernel(SGDGeneExpressionKernel, len(ps),
                            "p-gene-expression-kernel",
                            feat_to_i, self.src,
                            ["Gasch_2000_PMID_11102521",
                             "Spellman_1998_PMID_9843569"])

        # Compute the protein complex kernel
        self._cached_kernel(YeastProteinComplexKernel, len(ps),
                            "p-complex-kernel",
                            feat_to_i, self.src)

        # Compute the InterPro domain kernel
        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-kernel", ps, self.dst,
                            use_evalue = False)
        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-score-kernel", ps, self.dst,
                            use_evalue = True)

#        # Compute the profile kernel
#        self._cached_kernel(PSSMKernel, len(ps),
#                            "p-profile-kernel", p_to_i, self.dst,
#                            k = 4, threshold = 6.0)

    def _propagate_fun_on_dag(self, p_to_fun, dag):
        """WRITEME

        :param p_to_fun: map from protein IDs to GO term IDs.
        :param dag: the full GO DAG.
        :returns: map from protein IDs to ``GOTerm``'s.
        """
        propagated_p_to_fun = {}
        for p, term_ids in p_to_fun.iteritems():
            propagated_functions = set()
            for term_id in term_ids:
                if term_id == "": # This can occur due to empty annotations in SGD
                    continue
                paths = dag.paths_to_root(term_id)
                if paths == None:
                    print "Warning: protein '{}' is annotated with an undefined term '{}', skipping" \
                            .format(p, term_id)
                    continue
                for path in paths:
                    propagated_functions.update(path)
            propagated_p_to_fun[p] = propagated_functions
        return propagated_p_to_fun

    def _preprocess_dag(self, dag, ps, p_to_fun, aspects = None,
                        min_proteins_per_term = 30, max_depth = 5):
        """Processes the GO DAG by removing unwanted terms and adding bin terms.

        :param dag: the full GO DAG.
        :param ps: list of protein IDs.
        :param p_to_fun: map from protein IDs to ``GOTerm``'s.
        :param aspects: list of aspects to keep, in ``("BP", "CC", "MF")``.
        :param min_proteins_per_term: minimum number of proteins to keep a term.
        :param max_level: maximum term depth.
        :returns: the trimmed GO DAG, list of proteins, and protein-to-function map.
        """
        # TODO WRITEME
        return dag, ps, p_to_fun

    @staticmethod
    def _permute(l):
        pi = list(np.random.permutation(len(l)))
        return [l[pi[i]] for i in xrange(len(l))]

    def _get_folds(self, pps, p_to_terms, num_folds = 10):
        """Generates the folds.

        The interactions are split into folds according to their functions.
        No two protein-protein interactions can occur in two distinct folds.

        :param pps: list of interactions.
        :param p_to_terms: map between proteins and ``GOTerm``'s.
        :param num_folds: number of folds to generate (default: ``10``).
        :returns: a list of sets of protein pairs.
        """
        # Map from protein pairs to GO terms either protein is annotated with
        pp_to_terms = {(p1,p2): (set(p_to_terms[p1]) | set(p_to_terms[p2]))
                       for p1, p2 in pps}

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
        for term, pps in term_to_pps.iteritems():
            for p1, p2 in pps:
                assert (p2, p1) in pps

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

    def _write_sbr_data(self, ps, p_to_fun, dag, pp_pos, pp_neg):
        pass

    def _dump_dataset_statistics(self, ps, p_to_seq, p_to_fun, dag, prefix):
        """Dump a few dataset statistics."""

        # Draw the histogram of protein lengths
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = [len(p_to_seq[p]) for p in ps]
        ax.hist(x, range(min(x), max(x) + 50, 50), color = "red", alpha = 0.8)
        ax.yaxis.grid(True)
        fig.savefig(os.path.join(self.dst, "p-{}-length-hist.png".format(prefix)),
                    bbox_inches = "tight", pad_inches = 0)

        # Draw a histogram of the number of protein annotations per GO aspect
        # XXX we should find out proteins that have only the root annotated
        # XXX we should find out proteins with only *some* aspects

        # TODO: function co-occurrence
        # TODO: plot P(f|p)
        # TODO: plot P(p|f)
        # TODO: plot GO annotation depth
        # TODO: plot GO consistency w.r.t. one-path rule
        # TODO: plot GO consistency w.r.t. mutual exclusion
        # TODO: plot GO consistency w.r.t. interactions

    def run(self, min_sequence_len = 50, cdhit_threshold = 0.75):
        """Run the yeast prediction experiment."""

        # Retrieve the list of protein IDs (whose order will define the
        # structure of the kernels), and their mappings to feature IDs,
        # sequences and GO term annotations.
        ps          = self._cached(self._get_sgd_ids,
                                   "sgd_ids.pickle")
        p_to_feat   = self._cached(self._get_sgd_id_to_feat,
                                   "sgd_id_to_feat.pickle")
        p_to_seq    = self._cached(self._get_sgd_id_to_seq,
                                   "sgd_id_to_seq.pickle")
        p_to_fun    = self._cached(self._get_sgd_id_to_fun,
                                   "sgd_id_to_fun.pickle")
        print _cls(self), ": found {} proteins".format(len(ps))

        for p in ps:
            assert p in p_to_seq, "'{}' has no sequence".format(p)
            assert p in p_to_fun, "'{}' has no GO annotation".format(p)
            assert p in p_to_feat, "'{}' has no associated feature ID".format(p)

        # Load the GO OBO file
        dag = GODag(os.path.join(self.src, "GO", "go-basic.obo"))

        # TODO check if SGD GO annotations are GO annotations or GO slims, and
        # if the two map 1:1 to each other --- i.e. if we can use the GO IDs as
        # is.

        # Propagate the GO annotations to the root
        p_to_fun    = self._cached(self._propagate_fun_on_dag,
                                   "sgd_id_to_fun_propagated.pickle",
                                   p_to_fun, dag)

        # Dump the dataset statistics prior to any preprocessing
        self._dump_dataset_statistics(ps, p_to_seq, p_to_fun, dag, "raw")

        # Filter out sequences shorter than min_sequence_len residues, then
        # cluster the filtered sequences with CD-HIT
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= min_sequence_len, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), min_sequence_len)

        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, clusters = CDHit().run(filtered_p_seq, threshold = cdhit_threshold)
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(clusters), cdhit_threshold)

        filtered_ps = [list(cluster)[0][0] for cluster in clusters]

        # TODO preprocess the GO DAG and p_to_fun to remove (i) unwanted GO
        # aspects, (ii) terms with too few annotations, (iii) terms too deep in
        # the GO DAG, (iv) proteins with no annotations in the given aspects.

        # Query the high-quality protein-protein interactions, restricting
        # them to proteins in the filtered set
        pp_pos_hq, pp_pos_lq, pp_neg = self._get_sgd_pins(filtered_ps)

        # Process the GO DAG
        self._preprocess_dag(dag, filtered_ps, p_to_fun)

        # Dump the dataset statistics after the preprocessing
        self._dump_dataset_statistics(filtered_ps, p_to_seq, p_to_fun, dag,
                                      "preproc")

        # Create the folds
        print _cls(self), ": computing the folds"
        pp_folds = self._cached(self._get_folds, "folds",
                                pp_pos_hq | pp_neg, p_to_fun)

        # Compute the protein kernels
        print _cls(self), ": computing the protein kernels"
        self._compute_protein_kernels(filtered_ps, p_to_feat)

        # Write down the SBR input files
        print _cls(self), ": writing SBR input files"
        sbr_folds = self._cached(self._write_sbr_data, "sbr_folds",
                                 filtered_ps, p_to_fun, dag, pp_pos_hq, pp_neg)
