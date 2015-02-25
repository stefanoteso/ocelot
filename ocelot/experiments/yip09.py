# -*- coding: utf-8 -*-

from .base import _Experiment

import os
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

    def _get_sequences(self):
        """Reads the ORFs and their sequences from the endpoint."""
        ans = self.query("""
        SELECT ?p ?seq
        FROM <{default_graph}>
        WHERE {{
            ?p a ocelot:yip_protein .
            OPTIONAL {{
                ?feat a ocelot:sgd_feature .
                ?p owl:sameAs ?feat .
                ?orf a ocelot:sgd_id ;
                    ocelot:sgd_id_has_sequence ?seq .
                ?orf owl:sameAs ?feat .
            }}
        }}
        """)
        assert ans and len(ans) and "results" in ans
        p_to_seq = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p = bindings[u"p"].split(".")[-1]
            p_to_seq[p] = bindings[u"seq"]
        return p_to_seq

    def _get_ppis(self, symmetrize = False):
        """Reads the positive and negative PPIs from the endpoint.

        Please note that th Yip et al. dataset ignores the fact that the
        `interacts` predicate is symmetric --- read: their dataset is *not*
        symmetric.

        :param symmetrize: whether to symmetrize the interactions.
        """
        assert not symmetrize
        pos_ans = self.query("""
        SELECT ?p1 ?p2
        FROM <{default_graph}>
        WHERE {{
            ?p1 a ocelot:yip_protein .
            ?p2 a ocelot:yip_protein .
            ?p1 ocelot:yip_interacts_with ?p2 .
        }}
        """)
        neg_ans = self.query("""
        SELECT DISTINCT ?p1 ?p2
        FROM <{default_graph}>
        WHERE {{
            ?p1 a ocelot:yip_protein .
            ?p2 a ocelot:yip_protein .
            ?p1 ocelot:yip_not_interacts_with ?p2 .
        }}
        """)
        assert len(pos_ans["results"]["bindings"]) == 3201
        assert len(neg_ans["results"]["bindings"]) == 3201
        pos_ppi, neg_ppi = set(), set()
        for bindings in pos_ans["results"]["bindings"]:
            assert len(bindings) == 2
            pos_ppi.update([(self.cast(bindings[u"p1"]).split(".")[-1],
                             self.cast(bindings[u"p2"]).split(".")[-1])])
        for bindings in neg_ans["results"]["bindings"]:
            assert len(bindings) == 2
            neg_ppi.update([(self.cast(bindings[u"p1"]).split(".")[-1],
                             self.cast(bindings[u"p2"]).split(".")[-1])])
        assert len(pos_ppi) == 3201
        assert len(neg_ppi) == 3201
        assert len(pos_ppi & neg_ppi) == 0
        return pos_ppi, neg_ppi

    def _compute_ppi_y(self, pps):
        """Compute the y vector for protein-protein interactions.

        :param pps: sorted list of proteins.
        """
        print _cls(self), "retrieving protein-protein interactions"
        pos_ppi, neg_ppi = self._get_ppis()
        print " #ppi+={} #ppi-={}".format(len(pos_ppi), len(neg_ppi))

        path = os.path.join(self.dst, "ppi-y.txt")
        try:
            ys = np.loadtxt(path)
            assert len(ys) == (len(pos_ppi) + len(neg_ppi)), "y vector length mismatch"
        except Exception, e:
            print _cls(self), ":", e
            print _cls(self), "computing protein-protein y"
            ys = []
            for i, pp in enumerate(pps):
                if pp in pos_ppi:
                    ys.append(1.0)
                elif pp in neg_ppi:
                    ys.append(-1.0)
                else:
                    raise RuntimeError, "yip is angry with you '{}'".format(pp)
            ys = np.array(ys)
            np.savetxt(path, ys)
        return ys

    def _get_microarray_kernel(self, p_to_i, which = None):
        """Returns a kernel for gene expression.

        The data can be obtained `here`_.

        .. _here: ftp://downloads.yeastgenome.org/expression/microarray/all_spell_exprconn_pcl.tar.gz

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
        import numpy as np
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
        ans = self.query("""
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
        """)
        assert ans and len(ans) and "results" in ans
        context = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
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
        import numpy as np
        reader = PSSM()
        pssms = []
        for p, i in p_to_i.items():
            path = os.path.join(self.src, "yip09", "raw", "profiles",
                                "{}.ascii-pssm".format(p))
            print "loading {}".format(path)
            pssm = reader.read(path)
            pssms.append(pssm)
        return ProfileKernel(pssms)

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

        Note that the order in which objects are passed in is preserved when
        computing the kernels.

        XXX since all our kernel matrices are dense, this procedure is very
            memory intensive; we should definitely do something about that

        :param ps: list of protein identifiers
        :param pps: list of pairs of protein identifiers
        :param p_to_i: map between protein identifiers and indices
        """
        pp_indices = [ (p_to_i[p1], p_to_i[p2]) for p1, p2 in pps ]

        KERNEL_INFO = (
#            ("p-kernel-colocalization",
#                lambda: self._get_genetic_colocalization_kernel(p_to_i)),
            ("p-kernel-microarray",
                lambda: self._get_microarray_kernel(p_to_i,
                            which = ["Gasch_2000_PMID_11102521",
                                     "Spellman_1998_PMID_9843569"])),
#            ("p-kernel-microarray-all",
#                lambda: self._get_microarray_kerne(p_to_i)),
            ("p-kernel-complex",
                lambda: self._get_complex_kernel(p_to_i)),
            ("p-kernel-interpro-match-all",
                lambda: self._get_interpro_kernel(p_to_i, use_evalue = False)),
            ("p-kernel-interpro-weighted-all",
                lambda: self._get_interpro_kernel(p_to_i)),
            ("p-kernel-profile",
                lambda: self._get_profile_kernel(p_to_i)),
            ("p-kernel-yip",
                lambda: self._get_yip_kernels()[0]),
        )

        print _cls(self), ": computing protein kernels"
        kernels, pairwise_kernels = [], []
        for path, compute_kernel in KERNEL_INFO:
            path = os.path.join(self.dst, path)
            try:
                print _cls(self), ": loading '{}'".format(path)
                kernel = DummyKernel(path + ".txt", num = len(ps))
                assert kernel.is_psd(), "non-PSD kernel!"
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                kernel = compute_kernel()
                kernel.check_and_fixup(1e-10) 
                kernel.save(path + ".txt")
                kernel.draw(path + ".png")
            kernels.append(kernel)

            path += "-pairwise"
            try:
                print _cls(self), ": loading '{}".format(path)
                pairwise_kernel = DummyKernel(path + ".txt", num = len(pp_indices))
                assert kernel.is_psd(), "non-PSD kernel!"
            except Exception, e:
                print _cls(self), "|", e
                print _cls(self), ": computing '{}'".format(path)
                pairwise_kernel = PairwiseKernel(pp_indices, kernel)
                kernel.check_and_fixup(1e-10) 
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
        p_to_i = { p: i for i, p in enumerate(ps) }

        # The index of a pair of IDs is given by the order in which it appears
        # within the `goldPosProteinPairs.txt` etc. files, where all positives
        # are read first and all negatives follow. This is the same order in
        # which pairs appear in the pairwise kernels (XXX please double check)
        pp_to_i = { pp: i for i, pp in enumerate(pps) }

        # Retrieve the protein interactions
        pp_ys = self._compute_ppi_y(pps);

        # Compute the protein kernels
        p_kernels, pp_kernels = self._get_p_kernels(ps, pps, p_to_i)

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

