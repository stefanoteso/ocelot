# -*- coding: utf-8 -*-

from .base import _Experiment

import os
import numpy as np
import ocelot.ontology as O
from ocelot.services import _cls, iterate_csv
from ocelot.features import *
from ocelot.kernels import *
from ocelot.converters.yip09 import *

class YipExperiment(_Experiment):
    """Reproduce the experiment in [Yip09]_.

    We re-use the same labels and folds as the original dataset, but use the
    newly computed features.

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

    def _get_microarray_kernel(self, p_to_i, which = None):
        """Returns a kernel for gene expression.

        The data can be obtained `here <ftp://downloads.yeastgenome.org/expression/microarray/all_spell_exprconn_pcl.tar.gz>`_.
        """
        from glob import glob
        pcl = PCL()
        num = len(p_to_i)
        for dirname in glob(os.path.join(self.src, "yip09", "raw", "microarray", "*")):
            if not os.path.isdir(dirname):
                continue
            if which != None and not os.path.basename(dirname) in which:
                continue
            matrix_sum = np.zeros((num, num))
            for filename in glob(os.path.join(dirname, "*.pcl")):
                path = os.path.join(self.src, "yip09", "raw", "microarray",
                                    dirname, filename)
                print _cls(self), ": reading '{}'".format(path)
                orf_to_expression, num_conditions = pcl.read(path)
                levels = [ [0.0]*num_conditions for _ in xrange(len(p_to_i)) ]
                # XXX add support for bad/SGD IDs
                num_missing = 0
                for orf, index in p_to_i.items():
                    try:
                        levels[index] = orf_to_expression[orf]
                    except KeyError, key:
                        num_missing += 1
                if num_missing:
                    print "Warning: no expression level for '{}/{}' proteins".format(num_missing, len(p_to_i))
                matrix_sum += CorrelationKernel(levels).compute()
        return DummyKernel(matrix_sum)

    def _get_complex_kernel(self, p_to_i):
        """Computes diffusion kernels for raw protein complex data.

        The original experiment relies of Y2H [Ito00]_, [Uetz00]_ and TAP-MS
        [Gavin06]_, [Krogan06]_ datasets of protein complexes. Here we make
        use of the dataset described in [Pu08], which is supposed to be more
        complete and up-to-date.
        """
        import itertools as it
        # XXX use a weight matrix

        # Read the map between ORFs and complexes
        FIELDS = ("ORF", "_", "COMPLEX")
        complex_to_orfs = {}
        for row in iterate_csv(os.path.join(self.src, "yip09", "raw", "ppi", "CYC2008_complex.tab"),
                               delimiter = "\t", fieldnames = FIELDS):
            orf, cpx = row["ORF"], row["COMPLEX"]
            if not cpx in complex_to_orfs:
                complex_to_orfs[cpx] = set()
            complex_to_orfs[cpx].add(orf)

        # Now compute the adjacency matrix
        adj_matrix = np.zeros((len(p_to_i), len(p_to_i)))
        for _, orfs in complex_to_orfs.items():
            orf_indices = [p_to_i[orf] for orf in orfs if orf in p_to_i]
            for i, j in it.product(orf_indices, orf_indices):
                if i != j:
                    adj_matrix[i,j] = 1
        return DiffusionKernel(adj_matrix)

    def _get_interpro_kernel(self, p_to_i, allowed_sources = None,
                             use_evalue = False, default_score = 1.0):
        """Computes a set kernel over InterPro domain family hits.

        The InterPro files are read from:

        .. ``self.src/yip09/raw/interpro/${p}.iprscan5.tsv.txt``

        where ``${p}`` is the name of the protein.

        :param p_to_i: map protein names -> index in the kernel.
        :param allowed_sources: allowed InterPro sources, e.g. "Pfam"
        :param use_evalue: use negative logarithm of the hit e-value as
                           detection score.
        """
        interpro = InterProTSV()
        hits = [ None for _ in xrange(len(p_to_i)) ]
        for p, i in p_to_i.items():
            path = os.path.join(self.src, "yip09", "raw", "interpro",
                                "{}.iprscan5.tsv.txt".format(p))
            hit = interpro.read(path, allowed_sources)
            if use_evalue:
                # weight each hit by the negative log of its e-value
                for k, v in hit.items():
                    if v == None or v <= 0.0:
                        hit[k] = default_score
                    else:
                        hit[k] = -np.log(v)
            else:
                hit = set(hit.keys())
            hits[i] = hit
        if use_evalue:
            return SparseLinearKernel(hits)
        else:
            return SetKernel(hits)

    def _get_genetic_colocalization_kernel(self, p_to_i, gamma = 1e-9):
        """We use a simple Gaussian-of-differences kernel here.

        The distance between two proteins (or, rather, genes) is taken to
        be the distance between their centroids.

        Please note that we do distinguish between same-strand and
        different-strand proteins (i.e., their distances are computed the same
        way), while this may have a rather serious biological implications.

        The idea comes from [Lee03]_."""
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
        context = {}
        for bindings in self.iterquery(query, n = 5):
            p       = bindings[u"p"].split(".")[-1]
            chrom   = bindings[u"chrom"].split(".")[-1]
            strand  = bindings[u"strand"]
            start   = bindings[u"start"]
            stop    = bindings[u"stop"]
            assert strand in ("C", "W")
            if strand == "C":
                # XXX not really required, here just to not forget that 'W'
                # means increasing coords and 'C' means reverse.
                start, stop = stop, start
            context[p] = (chrom, min(start, stop) + 0.5 * np.fabs(start - stop))
        matrix = np.zeros((len(p_to_i), len(p_to_i)))
        for p, i in p_to_i.items():
            for q, j in p_to_i.items():
                if context[p][0] == context[q][0]:
                    matrix[i,j] = np.exp(-gamma * (context[p][1] - context[q][1])**2)
        return DummyKernel(matrix)

    def _get_profile_kernel(self, p_to_i):
        reader = PSSM(targets = ["residue", "nlog_condp"])
        pssms = []
        for p, i in p_to_i.items():
            path = os.path.join(self.src, "yip09", "raw", "profiles",
                                "{}.ascii-pssm".format(p))
            print "loading {}".format(path)
            pssm = reader.read(path)
            pssms.append(pssm)
        return ProfileKernel(pssms, k = 4, threshold = 6.0)

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

    def _get_p_kernels(self, ps, pps, p_to_i):
        """Computes all the kernels and pairwise kernels for proteins.

        The order in which objects are passed in is preserved.

        :param ps: list of protein identifiers
        :param pps: list of pairs of protein identifiers
        :param p_to_i: map from protein ID to protein index
        """
        INFOS = (
            ("p-kernel-colocalization",
                lambda: self._get_genetic_colocalization_kernel(p_to_i)),
            ("p-kernel-microarray",
                lambda: self._get_microarray_kernel(p_to_i,
                            which = ["Gasch_2000_PMID_11102521",
                                     "Spellman_1998_PMID_9843569"])),
            ("p-kernel-complex",
                lambda: self._get_complex_kernel(p_to_i)),
            ("p-kernel-interpro-match-all",
                lambda: self._get_interpro_kernel(p_to_i, use_evalue = False)),
            ("p-kernel-profile",
                lambda: self._get_profile_kernel(p_to_i)),
            ("p-kernel-yip",
                lambda: self._get_yip_kernels()[0]),
        )
        return self._get_x_kernels(INFOS, ps, pps, p_to_i, "protein")

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
        return self._get_x_kernels(INFOS, ds, dds, d_to_i, "domain")

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
        return self._get_x_kernels(INFOS, rs, rrs, r_to_i, "residue")

    def _get_x_kernels(self, infos, xs, xxs, x_to_i, name):
        """Helper for computing kernels and pairwise kernels."""
        print _cls(self), ": computing {} kernels".format(name)

        xx_indices = [(x_to_i[x1], x_to_i[x2]) for x1, x2 in xxs]

        kernels, pairwise_kernels = [], []
        for path, compute_kernel in infos:
            path = os.path.join(self.dst, path)
            try:
                print _cls(self), ": loading '{}'".format(path)
                kernel = DummyKernel(path + ".txt", num = len(xs))
                assert kernel.is_psd(), "non-PSD kernel!"
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                kernel = compute_kernel()
                if kernel == None:
                    continue
                kernel.check_and_fixup(1e-10) 
                kernel.save(path + ".txt")
                kernel.draw(path + ".png")
            if kernel:
                kernels.append(kernel)

            path += "-pairwise"
            try:
                print _cls(self), ": loading '{}".format(path)
                pairwise_kernel = DummyKernel(path + ".txt", num = len(xxs))
                assert kernel.is_psd(), "non-PSD kernel!"
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                pairwise_kernel = PairwiseKernel(xx_indices, kernel)
                assert not pairwise_kernel is None
                pairwise_kernel.check_and_fixup(1e-10) 
                pairwise_kernel.save(path + ".txt")
                pairwise_kernel.draw(path + ".png")
            pairwise_kernels.append(pairwise_kernel)

        return kernels, pairwise_kernels

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
        p_kernels, pp_kernels = self._get_p_kernels(ps, pps, p_to_i)
        d_kernels, dd_kernels = self._get_d_kernels(ds, dds, d_to_i, p_to_seq, d_to_pos, d_to_pfam)
        r_kernels, rr_kernels = self._get_r_kernels(rs, rrs, r_to_i, p_to_seq, r_to_pos)

        # Get the original fold splits and convert them to indices
        pp_ts_ids, dd_ts_ids, rr_ts_ids = converter.get_test_sets()

        pp_folds = []
        for ids in pp_ts_ids:
            ts_indices = [pp_to_i[(p1, p2)] for p1, p2 in ids]
            # XXX is this actually correct???
            tr_indices = sorted(set(range(len(pps))) - set(ts_indices))
            pp_folds.append((ts_indices, tr_indices))

        # Run the experiment
        print _cls(self), ": running"
        self._run_mkl(pp_ys, pp_kernels, pp_folds)

