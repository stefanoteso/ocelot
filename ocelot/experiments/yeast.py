# -*- coding: utf-8 -*-
from .base import _Experiment

import os
import numpy as np
from ocelot.services import _cls, CDHit
from ocelot.go import GODag, GOTerm

class YeastExperiment(_Experiment):
    """New experiment based on SGD and iPfam.

    This experiment is structured exactly the same as the Yip et al. [Yip09]_
    experiment, but performed on a different dataset based on SGD, BioGRID
    and iPfam.

    .. todo::

        Compute number of functions per protein and number of proteins per
        function; add to the GODag class.

    .. note::

        Somewhat unrelated. For up-to-date statistics on GO annotations,
        see `<http://geneontology.org/page/current-go-statistics>`_.

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
        return ids

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
        return sgd_id_to_feat

    def _get_sgd_pin(self, manual_only = False):
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
        return pp_pos

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

        query = """
        SELECT DISTINCT ?feat1 ?chain1 ?pfam1 ?feat2 ?chain2 ?pfam2 ?pdb
        FROM <{default_graph}>
        WHERE {{
            ?feat1 a ocelot:sgd_feature .
            ?hom1 a ocelot:sgd_pdb_homology ;
                ocelot:sgd_pdb_has_query ?feat1 ;
                ocelot:sgd_pdb_has_target ?chain1 ;
                ocelot:sgd_pdb_alignment ?evalue1 .
            ?region1 a ocelot:ipfam_region ;
                ocelot:ipfam_region_instance_of ?pfam1 ;
                ocelot:ipfam_region_occurs_in ?chain1 .

            ?feat2 a ocelot:sgd_feature .
            ?hom2 a ocelot:sgd_pdb_homology ;
                ocelot:sgd_pdb_has_query ?feat2 ;
                ocelot:sgd_pdb_has_target ?chain2 ;
                ocelot:sgd_pdb_alignment ?evalue2 .
            ?region2 a ocelot:ipfam_region ;
                ocelot:ipfam_region_instance_of ?pfam2 ;
                ocelot:ipfam_region_occurs_in ?chain2 .

            ?_ a ocelot:ipfam_region_int ;
                ocelot:ipfam_region_int_has_region ?region1 ;
                ocelot:ipfam_region_int_has_region ?region2 ;
                ocelot:ipfam_region_int_occurs_in ?pdb .
        }}"""
        dd_pos = set()
        for bindings in self.iterquery(query, n=7):
            print bindings
        sys.exit(1)
        return dd_pos

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

    def _propagate_fun_on_dag(self, p_to_fun, dag):
        propagated_p_to_fun = {}
        for p, functions in p_to_fun.iteritems():
            propagated_functions = set()
            for fun in functions:
                if fun == "": # This can occur due to empty annotations in SGD
                    continue
                paths = dag.paths_to_root(fun)
                if paths == None:
                    print "Warning: protein '{}' is annotated with an undefined term '{}', skipping" \
                            .format(p, fun)
                    continue
                for path in paths:
                    propagated_functions.update(path)
            propagated_p_to_fun[p] = propagated_functions
        return propagated_p_to_fun

    def _get_folds(ps, pos, num_folds = 10):
        raise NotImplementedError

    @staticmethod
    def _check_p_pp_are_sane(ps, pps):
        assert all((p, q) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"

    def run(self):
        print self._get_sgd_ipfam_din()
        """Run the yeast prediction experiment."""
        ps          = self._cached(self._get_sgd_ids,
                                   "sgd_ids.pickle")
        p_to_seq    = self._cached(self._get_sgd_id_to_seq,
                                   "sgd_id_to_seq.pickle")
        p_to_fun    = self._cached(self._get_sgd_id_to_fun,
                                   "sgd_id_to_fun.pickle")
        p_to_feat   = self._cached(self._get_sgd_id_to_feat,
                                   "sgd_id_to_feat.pickle")
        print _cls(self), ": found {} proteins".format(len(ps))

        # XXX curiously enough, some SGD proteins are annotated with GO terms
        # that are **not** part of goslim_yeast.obo... We use go-basic instead.
        dag = GODag(os.path.join(self.src, "GO", "go-basic.obo"))
        p_to_fun    = self._cached(self._propagate_fun_on_dag,
                                   "sgd_id_to_fun_propagated.pickle",
                                   p_to_fun, dag)

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
        pp_pos_hq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_hq.pickle",
                                 manual_only = True)
        density = float(len(pp_pos_hq)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} hi-quality PPIs (density = {})" \
                .format(len(pp_pos_hq), density)
        self._check_p_pp_are_sane(ps, pp_pos_hq)

        # Query the hq+lq protein-protein interactions
        pp_pos_lq = self._cached(self._get_sgd_pin,
                                 "sgd_id_interactions_lq.pickle",
                                 manual_only = False)
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

        # TODO retrieve dom-dom and res-res interaction instances

        # TODO create the folds, with redundancy reduced training sets
        pp_folds = self._cached(self._get_folds, "folds", ps, pp_pos_hq)
        print _cls(self), ": created the folds"

