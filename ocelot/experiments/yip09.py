# -*- coding: utf-8 -*-

from .base import _Experiment

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

    *References*

    .. [Yip09] Yip, Kim, McDermott, Gerstein, "Multi-level learning: improving
        the prediction of protein, domain and residue interactions by allowing
        information flow between levels", BMC Bioinformatics, 2009.
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

    def _get_microarray_kernel(self, p_to_i):
        """Returns a kernel for gene expression."""
        parts = {
            "Gasch_2000_PMID_11102521": (
                "2010.Gasch00_stationaryPhase(y12).flt.knn.avg.pcl",
                "2010.Gasch00_DTT(y14).flt.knn.avg.pcl",
                "2010.Gasch00_steadyState(y14).flt.knn.avg.pcl",
                "2010.Gasch00_steadyState(y13).flt.knn.avg.pcl",
                "2010.Gasch00_diamideTreatment.flt.knn.avg.pcl",
                "2010.Gasch00_HSto37.flt.knn.avg.pcl",
                "2010.Gasch00_HOtimeCourse.flt.knn.avg.pcl",
                "2010.Gasch00_Ndepletion.flt.knn.avg.pcl",
                "2010.Gasch00_menadione.flt.knn.avg.pcl",
                "2010.Gasch00_HSmild.flt.knn.avg.pcl",
                "2010.Gasch00_hypo-osmotic.flt.knn.avg.pcl",
                "2010.Gasch00_DTT(y13).flt.knn.avg.pcl",
                "2010.Gasch00_HS37-25.flt.knn.avg.pcl",
                "2010.Gasch00_hyper-osmotic.flt.knn.avg.pcl",
                "2010.Gasch00_HS25-37.flt.knn.avg.pcl",
                "2010.Gasch00_HS30-37.flt.knn.avg.pcl",
                "2010.Gasch00_adenineStarvation.flt.knn.avg.pcl",
                "2010.Gasch00_HS29-33.flt.knn.avg.pcl",
                "2010.Gasch00_carbonSources.flt.knn.avg.pcl",
                "2010.Gasch00_stationaryPhase(y14).flt.knn.avg.pcl",
            ),
            "Spellman_1998_PMID_9843569": (
                "2010.Spellman98_alphaFactor.flt.knn.avg.pcl",
                "2010.Spellman98_elutriation.flt.knn.avg.pcl",
                "2010.Spellman98_cdc15.flt.knn.avg.pcl",
            ),
        }
        pcl = PCL()
        num = len(p_to_i)
        for dirname in parts:
            matrix_sum = np.zeros((num, num))
            for filename in parts[dirname]:
                path = os.path.join(self.src, "yip09", "raw", "microarray",
                                    dirname, filename)
                print _cls(self), ": reading '{}'".format(path)
                orf_to_expression, num_conditions = PCL.read(path)
                levels = [ [0.0]*num_conditions for _ in xrange(len(p_to_i)) ]
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

        The Yip et al. dataset relies on a variety of sources, namely [Y2Ha]_,
        [Y2Hb]_, [TAPMSa]_, and [TAPMSb]_. However we ignore them and instead
        use [Cyc]_, which should be more complete and up-to-date.

        *References*

        .. [Y2Ha] Ito et al. "Toward a Protein-Protein Interaction Map of the
            Budding Yeast: A Comprehensive System to Examine Two-Hybrid
            Interactions in All Possible Combinations between the Yeast
            Proteins", PNAS, 2000.

        .. [Y2Hb] Uetz et al., "A Comprehensive Analysis of Protein-Protein
            Interactions in Saccharomyces cerevisiae", Nature, 2000.

        .. [TAPMSa] Gavin et al., "Proteome Survey Reveals Modularity of the
            Yeast Cell Machinery", Nature, 2006.

        .. [TAPMSb] Krogan et al, "Global Landscape of Protein Complexes in the
            Yeast Saccharomyces cerevisiae", Nature, 2006

        .. [Cyc] Pu et al., "Up-to-date catalogues of yeast protein complexes",
            NAR 2008
        """
        import itertools as it
        # XXX use a weight matrix

        # We ignore all of the following:
        #  - pnas_97_3_1143__5101Table2rev.xls
        #  - a wormhole
        #  - suppl[1234].xsl
        #  - 16554755id1088.tab.txt

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

    @staticmethod
    def _read_interpro(path, allowed_sources = None):
        """Reads an InterPro tabular-separated values file.

        :param path: path to the InterPro TSV file.
        :param allowed_sources: list of allowed domain sources (default: all).

        :returns: a list of ``(domain_id, evalue)`` pairs. ``evalue`` can be
                  None.
        """
        # XXX move to ocelot.services
        # XXX integrate iprscan5-urllib.py here; plumb `sequences` here
        FIELDS = ("QUERY_ID", "?2", "?3", "SOURCE_DB", "SOURCE_FAMILY",
                  "DESCRIPTION", "START", "STOP", "EVALUE", "?10", "DATE",
                  "IPR_FAMILY", "SHORT_DESCRIPTION", "GO_TERMS", "PATHWAYS")
        hits = {}
        for row in iterate_csv(path, delimiter="\t", fieldnames = FIELDS):
            if allowed_sources and not row["SOURCE_DB"] in allowed_sources:
                continue
            try:
                evalue = float(row["EVALUE"])
            except ValueError:
                evalue = None
            hits[row["SOURCE_DB"] + ":" + row["SOURCE_FAMILY"]] = evalue
        return hits

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
        hits = [ None for _ in xrange(len(p_to_i)) ]
        for p, i in p_to_i.items():
            path = os.path.join(self.src, "yip09", "raw", "interpro",
                                "{}.iprscan5.tsv.txt".format(p))
            hit = self._read_interpro(path, allowed_sources)
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

    def _compute_y_ppi(self, pos_ppi, neg_ppi, pps):
        print _cls(self), "computing protein-protein y"
        ys = []
        for i, pp in enumerate(pps):
            if pp in pos_ppi:
                ys.append(1.0)
            elif pp in neg_ppi:
                ys.append(-1.0)
            else:
                raise RuntimeError, "yip is angry with you '{}'".format(pp)
        return np.array(ys)

    def run(self):
        print _cls(self), ": retrieving entities"
        ps, ds, rs = Yip09Converter(self.src, None).get_entities()
        print " #ps={}, #ds={}, #rs={}".format(len(ps), len(ds), len(rs))

        print _cls(self), ": retrieving entity pairs"
        pps, dds, rrs = Yip09Converter(self.src, None).get_pairs()
        print " #pps={}, #dds={} #rrs={}".format(len(pps), len(dds), len(rrs))

        # This makes sure that proteins and protein pairs are sorted as in
        # the original Yip dataset.
        p_to_i = { p: i for i, p in enumerate(ps) }
        pp_indices = [ (p_to_i[p1], p_to_i[p2]) for p1, p2 in pps ]

        print _cls(self), ": retrieving protein sequences"
        p_to_seq = self._get_sequences()
        sequences = [p_to_seq[p] for p in ps]

        def linear_kernel_of(Features, sequences):
            phis = SequenceFeatures(Features).compute(sequences)
            return LinearKernel(phis)
        def set_kernel_of(Features, sequences):
            return None

        print _cls(self), ": computing kernels"
        KERNEL_INFO = (
            ("p-kernel-spectrum-k=2",
                lambda _: SpectrumKernel(_, kmin=2)),
            ("p-kernel-spectrum-k=3",
                lambda _: SpectrumKernel(_, kmin=3)),
            ("p-kernel-microarray",
                lambda _: self._get_microarray_kernel(p_to_i)),
            ("p-kernel-complex",
                lambda _: self._get_complex_kernel(p_to_i)),
            ("p-kernel-interpro-superfamily",
                lambda _: self._get_interpro_kernel(p_to_i,
                                                    allowed_sources = ("SUPERFAMILY",))),
            ("p-kernel-interpro-all",
                lambda _: self._get_interpro_kernel(p_to_i)),
        #    ("p-kernel-profile-it=2",
        #        lambda _: ProfileKernel(_, num_iterations=2)),
        )

        kernels, pairwise_kernels = [], []
        for path, func in KERNEL_INFO:
            path = os.path.join(self.dst, path)
            try:
                print _cls(self), ": loading '{}'".format(path)
                kernel = DummyKernel(path + ".txt")
                matrix = kernel.compute()
                assert matrix.shape == (len(ps), len(ps))
            except Exception, e:
                print _cls(self), ":", e
                print _cls(self), ": computing '{}'".format(path)
                kernel = func(sequences)
                if kernel != None:
                    matrix = kernel.compute()
                    assert (matrix.T == matrix).all(), "not symmetric!"
                    ls = np.linalg.eigvalsh(matrix)
                    if ls[0] < 0.0:
                        assert ls[0] > -1e-10, "matrix is too non-PSD"
                        print _cls(self), ": preconditioning by 10*(ls[0] = '{}')".format(ls[0])
                        matrix += np.identity(matrix.shape[0]) * -10.0 * ls[0]
                        assert np.linalg.eigvalsh(matrix)[0] >= 0, "the gods are playful today"
                    print _cls(self), ": saving '{}'".format(path)
                    np.savetxt(path + ".txt", matrix)
                    kernel.draw(path + ".png")
            if kernel:
                del kernel

            #try:
            #    print _cls(self), ": loading pairwise '{}'".format(path)
            #    pairwise_kernel = DummyKernel(path + "-pairwise.txt")
            #    pairwise_matrix = pairwise_kernel.compute()
            #    assert pairwise_matrix.shape == (len(pps), len(pps))
            #except:
            #    print _cls(self), ": computing pairwise '{}'".format(path)
            #    pairwise_kernel = PairwiseKernel(pp_indices, kernel)
            #    pairwise_matrix = pairwise_kernel.compute()
            #    np.savetxt(path + "-pairwise.txt", pairwise_matrix)
            #    pairwise_kernel.draw(path + "-pairwise.png")
            #if pairwise_kernel:
            #    del pairwise_kernel

        print _cls(self), ": retrieving protein interactions"
        pos_ppi, neg_ppi = self._get_ppis()
        print " #ppi+={} #ppi-={}".format(len(pos_ppi), len(neg_ppi))

        print _cls(self), ": computing labels for ppi"
        ys = self._compute_y_ppi(pos_ppi, neg_ppi, pps)
        np.savetxt(os.path.join(self.dst, "ppi-y.txt"), ys)