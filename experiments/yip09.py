# -*- coding: utf-8 -*-

import os
import numpy as np
from ocelot.services import _cls
from ocelot.kernels import *
from ocelot.converters.yip09 import *
from ocelot import Stage, Experiment
from .yeast import *

_KERNEL_RELPATHS_PAIRWISE = (
    "p-colocalization-kernel-pairwise",
    "p-gene-expression-kernel-pairwise",
    "p-complex-kernel-pairwise",
    "p-interpro-kernel-pairwise",
    "p-interpro-score-kernel-pairwise",
    "p-profile-kernel-pairwise",
    "p-random-kernel-pairwise",
)

class YipExperiment(Experiment):
    """Reproduce the experiment in [Yip09]_.

    We re-use the same labels and folds as the original dataset, but use a new
    set of kernels (computed along the same lines as the in the original
    dataset, plus some).

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    def __init__(self, *args, **kwargs):
        super(YipExperiment, self).__init__(*args, **kwargs)

    def _get_p_to_seq(self):
        query = """
        SELECT ?yip_p ?seq
        FROM <{default_graph}>
        WHERE {{
            ?yip_p a ocelot:yip_protein ;
                owl:sameAs ?sgd_feat .
            ?sgd_feat a ocelot:sgd_feature .
            ?sgd_id a ocelot:sgd_id ;
                owl:sameAs ?sgd_feat ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_sequence ?seq .
        }}
        """
        p_to_seq = {binding[u"yip_p"].split(".")[-1]: binding[u"seq"]
                    for binding in self.iterquery(query, n = 2)}
        return p_to_seq

    def _get_p_to_context(self):
        query = """
        SELECT ?p ?chrom ?strand ?start ?stop
        FROM <{default_graph}>
        WHERE {{
            ?p a ocelot:yip_protein ;
                owl:sameAs ?feat .
            ?feat a ocelot:sgd_feature .
            ?id a ocelot:sgd_id ;
                owl:sameAs ?feat .
            ?id ocelot:sgd_id_in_chromosome ?chrom .
            ?id ocelot:sgd_id_in_strand ?strand .
            ?id ocelot:sgd_id_starts_at ?start .
            ?id ocelot:sgd_id_stops_at ?stop .
        }}
        """
        p_to_context = {}
        for bindings in self.iterquery(query, n = 5):
            p       = bindings[u"p"].split(".")[-1]
            chrom   = bindings[u"chrom"].split(".")[-1]
            strand  = bindings[u"strand"]
            start   = bindings[u"start"]
            stop    = bindings[u"stop"]
            assert strand in ("C", "W")
            p_to_context[p] = \
                (chrom, min(start, stop) + 0.5 * np.fabs(start - stop))
        return p_to_context

    def _get_d_r_info(self):
        # Note that all positions start from 1.
        # XXX this only fetches domains that *do* have child residues in the
        # yip dataset; one exception is YAL038W, which has no residues -- and
        # is filtered out by the query. Will fix it later.
        query = """
        SELECT ?yip_d ?yip_r ?d_pos0 ?d_pos1 ?r_pos ?pfam
        FROM <{default_graph}>
        WHERE {{
            # a yip protein
            ?yip_p a ocelot:yip_protein ;
                ocelot:yip_parent_of ?yip_d ;
                owl:sameAs ?sgd_feat .

            # a yip domain, child of the yip protein
            ?yip_d a ocelot:yip_domain ;
                ocelot:yip_parent_of ?yip_r ;
                ocelot:yip_instance_of ?pfam .

            # a yip residue, child of the yip residue
            ?yip_r a ocelot:yip_residue ;
                ocelot:yip_residue_has_position ?r_pos .

            # map the yip protein to an SGD ORF
            ?sgd_feat a ocelot:sgd_feature ;
                ocelot:sgd_feature_has_ipr_hit ?ipr_hit .
            ?sgd_id a ocelot:sgd_id ;
                owl:sameAs ?sgd_feat ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF .

            # get the InterPro hits for the SGD ORF that match the protein
            # family of the yip domain
            ?ipr_hit a ocelot:sgd_ipr_hit ;
                ocelot:sgd_ipr_hit_is_true true ;
                ocelot:sgd_ipr_hit_has_method "Pfam" ;
                ocelot:sgd_ipr_hit_has_db_id ?pfam ;
                ocelot:sgd_ipr_hit_starts_at ?d_pos0 ;
                ocelot:sgd_ipr_hit_stops_at ?d_pos1 .

            FILTER((?r_pos >= ?d_pos0) && (?r_pos <= ?d_pos1))
        }}
        """
        # The yip dataset does not contain duplicates; that is, the child
        # domains of a protein always have *distinct* Pfam IDs. However, we
        # should make sure that the domains are also unique.
        d_to_pos, d_to_pfam, r_to_pos = {}, {}, {}
        for binding in self.iterquery(query):
            d = binding[u"yip_d"].split(".")[-1]
            r = binding[u"yip_r"].split(".")[-1]
            d_to_pos[d] = (binding[u"d_pos0"], binding[u"d_pos1"])
            r_to_pos[r] = binding[u"r_pos"]
            d_to_pfam[d] = binding[u"pfam"]
        return d_to_pos, r_to_pos, d_to_pfam

    def _get_ppis(self, symmetrize = False):
        """Reads the positive and negative PPIs from the endpoint.

        Please note that th Yip et al. dataset ignores the fact that the
        `interacts` predicate is symmetric --- read: their dataset is *not*
        symmetric.

        :param symmetrize: whether to symmetrize the interactions.
        """
        assert not symmetrize
        pos_query = """
        SELECT ?p1 ?p2
        FROM <{default_graph}>
        WHERE {{
            ?p1 a ocelot:yip_protein .
            ?p2 a ocelot:yip_protein .
            ?p1 ocelot:yip_interacts_with ?p2 .
        }}
        """
        neg_query = """
        SELECT DISTINCT ?p1 ?p2
        FROM <{default_graph}>
        WHERE {{
            ?p1 a ocelot:yip_protein .
            ?p2 a ocelot:yip_protein .
            ?p1 ocelot:yip_not_interacts_with ?p2 .
        }}
        """
        def binding_to_pair(binding):
            return binding[u"p1"].split(".")[-1], binding[u"p2"].split(".")[-1]
        pos_ppi = set(binding_to_pair(binding)
                      for binding in self.iterquery(pos_query, n = 2))
        neg_ppi = set(binding_to_pair(binding)
                      for binding in self.iterquery(neg_query, n = 2))
        assert len(pos_ppi) == 3201
        assert len(neg_ppi) == 3201
        assert len(pos_ppi & neg_ppi) == 0
        return pos_ppi, neg_ppi

    def _get_ddis(self, symmetrize = False):
        """Reads the positive and negative DDIs from the endpoint."""
        assert not symmetrize
        pos_query = """
        SELECT ?d1 ?d2
        FROM <{default_graph}>
        WHERE {{
            ?d1 a ocelot:yip_domain .
            ?d2 a ocelot:yip_domain .
            ?d1 ocelot:yip_interacts_with ?d2 .
        }}
        """
        neg_query = """
        SELECT ?d1 ?d2
        FROM <{default_graph}>
        WHERE {{
            ?d1 a ocelot:yip_domain .
            ?d2 a ocelot:yip_domain .
            ?d1 ocelot:yip_not_interacts_with ?d2 .
        }}
        """
        def binding_to_pair(binding):
            return tuple(binding[u"d1"].split(".")[-1].split("_")), \
                   tuple(binding[u"d2"].split(".")[-1].split("_"))
        pos_ddi = set(binding_to_pair(binding)
                      for binding in self.iterquery(pos_query, n = 2))
        neg_ddi = set(binding_to_pair(binding)
                      for binding in self.iterquery(neg_query, n = 2))
        assert len(pos_ddi) == 422, len(pos_ddi)
        assert len(neg_ddi) == 422, len(neg_ddi)
        assert len(pos_ddi & neg_ddi) == 0
        return pos_ddi, neg_ddi

    def _get_rris(self, symmetrize = False):
        """Reads the positive and negative RRIs from the endpoint."""
        pos_query = """
        SELECT ?r1 ?r2
        FROM <{default_graph}>
        WHERE {{
            ?r1 a ocelot:yip_residue .
            ?r2 a ocelot:yip_residue .
            ?r1 ocelot:yip_interacts_with ?r2 .
        }}
        """
        neg_query = """
        SELECT ?r1 ?r2
        FROM <{default_graph}>
        WHERE {{
            ?r1 a ocelot:yip_residue .
            ?r2 a ocelot:yip_residue .
            ?r1 ocelot:yip_not_interacts_with ?r2 .
        }}
        """
        def binding_to_pair(binding):
            return tuple(binding[u"r1"].split(".")[-1].split("_")), \
                   tuple(binding[u"r2"].split(".")[-1].split("_"))
        pos_ddi = set(binding_to_pair(binding)
                      for binding in self.iterquery(pos_query, n = 2))
        neg_ddi = set(binding_to_pair(binding)
                      for binding in self.iterquery(neg_query, n = 2))
        assert len(pos_ddi) == 2000, len(pos_ddi)
        assert len(neg_ddi) == 2000, len(neg_ddi)
        assert len(pos_ddi & neg_ddi) == 0
        return pos_ddi, neg_ddi

    def _compute_ppi_y(self, pps):
        return self._compute_xxi_y(pps, self._get_ppis, "ppi")

    def _compute_ddi_y(self, dds):
        return self._compute_xxi_y(dds, self._get_ddis, "ddi")

    def _compute_rri_y(self, rrs):
        return self._compute_xxi_y(rrs, self._get_rris, "rri")

    def _compute_xxi_y(self, pairs, get_interactions, name):
        print _cls(self), ": retrieving {} interactions".format(name)
        pos_ints, neg_ints = get_interactions()
        print " #{}+={} #{}-={}".format(name, len(pos_ints), name, len(neg_ints))

        path = os.path.join(self.dst, "{}-y.txt".format(name))
        try:
            ys = np.loadtxt(path)
            assert len(ys) == (len(pos_ints) + len(neg_ints))
        except Exception, e:
            print _cls(self), ":", e
            print _cls(self), ": computing {} y".format(name)
            ys = []
            for i, pair in enumerate(pairs):
                if pair in pos_ints:
                    ys.append(1.0)
                elif pair in neg_ints:
                    ys.append(-1.0)
                else:
                    raise ValueError, "yip is angry with you '{}'".format(pair)
            ys = np.array(ys)
            np.savetxt(path, ys)
        return ys

    def _get_domain_ptcorr_kernel(self):
        # TODO: for each pfam pair, compute the phylogenetic tree correlation,
        #       giving a real matrix; empirical kernel map
        pass

    def _get_domain_family_frequency_kernel(self):
        # TODO: for each pfam, count the number of proteins that have a domain
        #       in that pfam; poly 3
        # TODO: for each pfam, count the number of proteins that have only
        #       domains of that pfam; poly 3
        pass

    def _get_domain_frequency_kernel(self):
        # TODO: for each domain, count the number of domains in the parent
        #       protein; poly 3
        pass

    def _get_residue_profile_kernel(self):
        pass

    def _get_residue_ss_kernel(self):
        pass

    def _get_residue_sa_kernel(self):
        pass

    def _get_yip_kernels(self):
        try:
            return self.yip_p_kernel, self.yip_d_kernel, self.yip_r_kernel
        except:
            converter = Yip09Converter(self.src, None, basename = "yip09")
            self.yip_p_kernel, self.yip_d_kernel, self.yip_r_kernel = \
                converter.get_kernels()
        return self.yip_p_kernel, self.yip_d_kernel, self.yip_r_kernel

    def _get_p_kernels(self, ps, pps, p_to_i, p_to_seq):
        """Computes all the kernels and pairwise kernels for proteins.

        The order in which objects are passed in is preserved.

        :param ps: list of protein IDs.
        :param pps: list of pairs of protein IDs.
        :param p_to_i: map from protein ID to protein index.
        :param p_to_seq: map from protein ID to protein sequence.
        """
        p_to_i = {p: i for i, p in enumerate(ps)}
        pp_indices = [(p_to_i[p1], p_to_i[p2]) for p1, p2 in pps]

        p_to_context = self._cached(self._get_p_to_context,
                                    "p_to_context")
        contexts = [p_to_context[p] for p in ps]
        self._cached_kernel(ColocalizationKernel, len(ps),
                            "p-colocalization-kernel",
                            contexts, gamma = 1.0,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(SGDGeneExpressionKernel, len(ps),
                            "p-gene-expression-kernel",
                            p_to_i, self.src,
                            ["Gasch_2000_PMID_11102521",
                             "Spellman_1998_PMID_9843569"],
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(YeastProteinComplexKernel, len(ps),
                            "p-complex-kernel",
                            p_to_i, self.src,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-kernel",
                            ps, self.dst, use_evalue = False,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(InterProKernel, len(ps),
                            "p-interpro-score-kernel",
                            ps, self.dst, use_evalue = True,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(PSSMKernel, len(ps),
                            "p-profile-kernel",
                            ps, p_to_seq, self.dst, k = 4, threshold = 6.0,
                            do_pairwise = True, pairs = pp_indices)

        self._cached_kernel(RandomKernel, len(ps),
                            "p-random-kernel",
                            ps, self.dst,
                            do_pairwise = True, pairs = pp_indices)

        print _cls(self), ": computing average kernel"
        avg_matrix = self._compute_average_kernel(_KERNEL_RELPATHS_PAIRWISE[:-1])
        kernel = DummyKernel(avg_matrix)
        kernel.check_and_fixup(1e-10)
        kernel.save(os.path.join(self.dst, "p-average-kernel-pairwise.txt"))
        kernel.draw(os.path.join(self.dst, "p-average-kernel-pairwise.png"))

    def _get_d_kernels(self, ds, dds, d_to_i, d_to_pos, d_to_pfam):
        """Computes all the kernels and pairwise kernels for domains.

        The order in which objects are passed in is preserved.

        .. todo::

            We are missing two domain kernels still. They need per-fold
            computation though, so they are not trivial to stitch in.

        :param ds: list of domain IDs.
        :param dds: list of domain ID pairs.
        :param d_to_i: map from domain ID to kernel index.
        :param p_to_seq: map from protein ID to protein sequence.
        :param d_to_pos: map from domain ID to position within the protein.
        """
        INFOS = (
            ("d-kernel-ptcorr",
                lambda: self._get_domain_ptcorr_kernel()),
            ("d-kernel-family-freq",
                lambda: self._get_domain_family_freq_kernel()),
            ("d-kernel-parent-freq",
                lambda: self._get_domain_parent_freq_kernel()),
        )
        return self._compute_kernels(INFOS, ds, dds)

    def _get_r_kernels(self, rs, rrs, r_to_i, p_to_seq, r_to_pos):
        """Computes all the kernels and pairwise kernels for residues.

        The order in which objects are passed in is preserved.

        :param rs: list of residue IDs.
        :param rrs: list of residue ID pairs.
        :param r_to_i: map from residue ID to kernel index.
        :param p_to_seq: map from protein ID to protein sequence.
        :param r_to_pos: map from residue ID to position within the protein.
        """
        INFOS = (
            ("r-kernel-profile",
                lambda: self._get_residue_profile_kernel()),
            ("r-kernel-ss",
                lambda: self._get_residue_ss_kernel()),
            ("r-kernel-sa",
                lambda: self._get_residue_sa_kernel()),
        )
        return self._compute_kernels(INFOS, rs, rrs)

    def _test_pp_kernels(self, folds, ys):
        from pprint import pprint
        CS = 10.0**np.linspace(-4, 4, 3)
        for relpath in _KERNEL_RELPATHS_PAIRWISE:
            gram = self._load_kernel(relpath)
            for c in CS:
                results = self.eval_svm(folds, ys, gram, c = c)
                for sup, f1, pr, rc, auc in results:
                    print "{} -> {} : F1={} Pr={} Rc={} AUC={}".format(c, sup, f1, pr, rc, auc)

    def run(self):
        """Run the Yip et al. experiment replica."""

        # Not everything is converted to RDF; namely, the *order* in which
        # proteins, domains, residues (and their interactions) should appear
        # is not. Here we read those from the original dataset itself.
        # XXX this is only required because we may want to use our kernels
        # as drop-in replacements in the existing SBR/yip experiment.
        converter = Yip09Converter(self.src, None, basename = "yip09")

        print _cls(self), ": retrieving entities and entity pairs"
        ps, ds, rs = converter.get_entities()
        pps, dds, rrs = converter.get_pairs()
        print " #ps={}, #ds={}, #rs={}".format(len(ps), len(ds), len(rs))
        print " #pps={}, #dds={} #rrs={}".format(len(pps), len(dds), len(rrs))

        # The index of an ID is given by the order in which it appears within
        # the `proteins.txt`, `domains.txt` or `residues.txt` files. This is
        # the same order in which entities appear in the non-pairwise kernels.
        # The index of a pair of IDs is given by the order in which it appears
        # within the `goldPosProteinPairs.txt` etc. files, where all positives
        # are read first and all negatives follow. This is the same order in
        # which pairs appear in the pairwise kernels (XXX please double check)
        def f(coll):
            return {x: i for i, x in enumerate(coll)}

        p_to_i, pp_to_i = f(ps), f(pps)
        d_to_i, dd_to_i = f(ds), f(dds)
        r_to_i, rr_to_i = f(rs), f(rrs)

        # Retrieve the interactions
        pp_ys = self._compute_ppi_y(pps)
        dd_ys = self._compute_ddi_y(dds)
        rr_ys = self._compute_rri_y(rrs)

        # Retrieve the protein sequences
        p_to_seq = self._get_p_to_seq()

        # Compute the domain/residue positions (in a best effort fashion)
        d_to_pos, r_to_pos, d_to_pfam = self._get_d_r_info()

        # Compute the protein kernels
        self._get_p_kernels(ps, pps, p_to_i, p_to_seq)
        self._get_d_kernels(ds, dds, d_to_i, p_to_seq, d_to_pos, d_to_pfam)
        self._get_r_kernels(rs, rrs, r_to_i, p_to_seq, r_to_pos)

        # Get the original fold splits and convert them to indices
        pp_folds, dd_folds, rr_folds = converter.get_test_sets()

        pp_folds_i = []
        for k in xrange(len(pp_folds)):
            pp_test_indices = [pp_to_i[(p1,p2)] for p1,p2 in pp_folds[k]]
            pp_train_indices = []
            for l in xrange(len(pp_folds)):
                pp_train_indices.extend([pp_to_i[(p1,p2)] for p1,p2 in pp_folds[k]
                                         if k != l])
            pp_folds_i.append((pp_test_indices, pp_train_indices))

        # Test each kernel independently using an SVM
        self._test_pp_kernels(pp_folds_i, pp_ys)
