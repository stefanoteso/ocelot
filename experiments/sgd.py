# -*- coding: utf-8 -*-

from . import _Experiment

import os
import numpy as np
import itertools as it
from glob import glob
from ocelot.services import _cls, iterate_csv, CDHit, PCL
from ocelot.go import GODag, GOTerm
from ocelot.kernels import *

class SGDGeneExpressionKernel(Kernel):
    """A yeast-specific correlation kernel on gene expression data.

    It uses the microarray files located at::

        $src/yeast/microarray/*.pcl

    :param p_to_i: map between ORF feature names to indices in the kernel matrix.
    :param src: path to the data source directory.
    :param experiments: list of microarray experiments to use.
    """
    def __init__(self, p_to_i, src, experiments, *args, **kwargs):
        self._wc = os.path.join(src, "SGD", "microarray", "*", "*.pcl")
        self._experiments = experiments
        super(SGDGeneExpressionKernel, self).__init__(p_to_i, *args, **kwargs)

    def _get_pcl_paths(self):
        for path in glob(self._wc):
            if not any(exp in path for exp in self._experiments):
                continue
            yield path

    def _compute_all(self):
        pcl = PCL()
        p_to_i = self._entities
        matrix = np.zeros((len(self), len(self)))
        for path in self._get_pcl_paths():
            print _cls(self), ": processing '{}'".format(path)
            p_to_levels, num_conditions = pcl.read(path)
            exp_levels, num_missing = [], 0
            for p, i in sorted(p_to_i.iteritems(), key = lambda p_i: p_i[1]):
                try:
                    p_levels = p_to_levels[p]
                except:
                    p_levels = np.zeros((num_conditions,))
                    num_missing += 1
                exp_levels.append(p_levels)
            if num_missing > 0:
                 print _cls(self), ": '{}' has no measurements for '{}/{}' proteins" \
                                    .format(path, num_missing, len(self))
            matrix += CorrelationKernel(exp_levels, do_normalize = self._do_normalize).compute()
            break
        return matrix

class YeastProteinComplexKernel(Kernel):
    """A yeast-specific diffusion kernel on protein complexes.

    It uses the protein complex file located at::

        $src/yeast/ppi/CYC2008_complex.tab

    :param p_to_i: map between ORF feature names to indices in the kernel matrix.
    :param src: path to the data source directory.
    :param gamma: parameter of the diffusion kernel.
    """
    def __init__(self, p_to_i, src, *args, **kwargs):
        self._path = os.path.join(src, "yeast", "ppi", "CYC2008_complex.tab")
        super(YeastProteinComplexKernel, self).__init__(p_to_i, *args, **kwargs)

    def _read_complex_to_orf(self):
        FIELDS = ("ORF", "_", "COMPLEX")
        complex_to_orfs = {}
        for row in iterate_csv(self._path, delimiter = "\t", num_skip = 1,
                               fieldnames = FIELDS):
            orf, cplex = row["ORF"], row["COMPLEX"]
            if not cplex in complex_to_orfs:
                complex_to_orfs[cplex] = set()
            complex_to_orfs[cplex].add(orf)
        return complex_to_orfs

    def _compute_all(self):
        p_to_i = self._entities

        complex_to_orfs = self._read_complex_to_orf()

        adjmatrix = np.zeros((len(self), len(self)))
        for _, orfs in complex_to_orfs.iteritems():
            orf_indices = [p_to_i[orf] for orf in orfs if orf in p_to_i]
            for i, j in it.product(orf_indices, orf_indices):
                if i != j:
                    adjmatrix[i,j] = 1
        return DiffusionKernel(adjmatrix).compute()

class InterProKernel(Kernel):
    """A simple domain kernel built around InterPro.

    :param p_to_i: map between protein identifiers and kernel index.
    :param cache_path: path to the directory holding the interpro output.
    :param allowed_sources: list of allowed domain providers (default: None).
    :param use_evalue: whether to use the evalue to score the predictions.
    """
    def __init__(self, p_to_i, cache_path, allowed_sources = None,
                 use_evalue = False, *args, **kwargs):
        self._cache_path = cache_path
        self._allowed_providers = allowed_providers
        self._use_evalue = use_evalue
        super(InterProKernel, self).__init__(p_to_i, *args, **kwargs)

    def _path(self, p):
        return os.path.join(self.cache_path, "interpro", "{}.f.tsv.txt".format(p))

    def _retrieve_interpro_predictions(self):
        raise NotImplementedError

    def _compute_all(self):
        self._retrieve_interpro_predictions()

        parser = InterProTSV()

        all_hits = []
        for p, i in sorted(p_to_i.items(), key = lambda p_i: p_i[1]):
            domain_to_evalue = interpro.read(self._path(p), allowed_sources)
            if not self._use_evalue:
                hits = set(domain_to_evalue.keys())
            else:
                hits = {domain: (-np.log(evalue) if evalue is None else DEFAULT_SCORE)
                        for domain, evalue in domain_to_evalue.iteritems()}
            all_hits.append(hits)

        if not self._use_evalue:
            return SetKernel(all_hits).compute()
        else:
            return SparseLinearKernel(all_hits).compute()

class PSSMKernel(ProfileKernel):
    """A simple wrapper around the profile kernel for strings [Kuang04]_.

    It takes care of generating the PSSM profiles.

    :param p_to_i: WRITEME
    :param cache_path: WRITEME
    :param num_iterations: WRITEME (default: 2)

    All remaining options are passed to the underlying ``ProfileKernel``.
    """
    def __init__(self, p_to_i, cache_path, num_iterations = 2, *args, **kwargs):
        self._cache_path = cache_path
        super(PSSMKernel, self).__init__(p_to_i, *args, **kwargs)

    def _get_pssm_path(self, p):
        return os.path.join(self.cache_path, "pssm", "{}.ascii-pssm".format(p))

    def _compute_pssms(self):
        raise NotImplementedError

    def _compute_all(self):
        self._compute_pssms()

        reader = PSSM(targets = ("residue", "nlog_condp"))
        pssms = []
        for p, i in sorted(p_to_i.items(), key = lambda p_i: p_i[1]):
            pssms.append(reader.read(self._get_pssm_path(p)))
        self._entities = pssms
        return super(PSSMKernel, self)._compute_all()

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

#        # Compute the InterPro domain kernel
#        self._cached_kernel(InterProKernel, len(ps),
#                            "p-interpro-kernel", p_to_i, self.dst)
#
#        # Compute the profile kernel
#        self._cached_kernel(PSSMKernel, len(ps),
#                            "p-profile-kernel", p_to_i, self.dst,
#                            k = 4, threshold = 6.0)

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

    def _draw_dataset_statistics(self, ps, p_to_seq, p_to_fun):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = [len(p_to_seq[p]) for p in ps]
        ax.hist(x, range(min(x), max(x) + 50, 50), color = "red", alpha = 0.8)
        ax.yaxis.grid(True)
        fig.savefig(os.path.join(self.dst, "p-length-hist.png"),
                    bbox_inches = "tight", pad_inches = 0)

        # TODO: function co-occurrence
        # TODO: plot P(f|p)
        # TODO: plot P(p|f)
        # TODO: plot GO annotation depth
        # TODO: plot GO consistency w.r.t. one-path rule
        # TODO: plot GO consistency w.r.t. mutual exclusion
        # TODO: plot GO consistency w.r.t. interactions

    @staticmethod
    def _check_p_pp_are_sane(ps, pps):
        assert all((p, q) in pps for (p, q) in pps), \
            "pairs are not symmetric"
        assert all(p in ps for p, _ in pps), \
            "singletons and pairs do not match"

    def run(self, min_sequence_len = 50, cdhit_threshold = 0.75):
        """Run the yeast prediction experiment."""

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

        # XXX curiously enough, some SGD proteins are annotated with GO terms
        # that are **not** part of goslim_yeast.obo... We use go-basic instead.
        dag = GODag(os.path.join(self.src, "GO", "go-basic.obo"))
        p_to_fun    = self._cached(self._propagate_fun_on_dag,
                                   "sgd_id_to_fun_propagated.pickle",
                                   p_to_fun, dag)

        self._draw_dataset_statistics(ps, p_to_seq, p_to_fun)

        # Filter out sequences shorter than min_sequence_len residues
        filtered_ps = filter(lambda p: len(p_to_seq[p]) >= min_sequence_len, ps)
        print _cls(self), ": found {} proteins with at least {} residues" \
                .format(len(ps), min_sequence_len)

        # Cluster sequences with CD-HIT
        filtered_p_seq = zip(filtered_ps, [p_to_seq[p] for p in filtered_ps])
        _, clusters = CDHit().run(filtered_p_seq, threshold = cdhit_threshold)
        print _cls(self), ": found {} clusters of proteins at CD-HIT threshold {}" \
                .format(len(clusters), cdhit_threshold)

        filtered_ps = [list(cluster)[0][0] for cluster in clusters]

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

        # Compute the protein kernels
        self._compute_protein_kernels(filtered_ps, p_to_feat)

        # TODO retrieve dom-dom and res-res interaction instances

        # TODO create the folds, with redundancy reduced training sets
        pp_folds = self._cached(self._get_folds, "folds", ps, pp_pos_hq)
        print _cls(self), ": created the folds"


