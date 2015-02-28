# -*- coding: utf-8 -*-

from .base import _Experiment

from ocelot.services import _cls

class YeastExperiment(_Experiment):
    """New experiment based on SGD and iPfam.

    This experiment is structured exactly the same as the Yip et al. [Yip09]_
    experiment, except on a new dataset.

    :param src: source directory of all databases.
    :param dst: destination directory for the results.
    :param subtargets: prediction targets.
    :param endpoint: URI of the SPARQL endpoint.
    :param default_graph: URI of the default graph.
    """
    def __init__(self, *args, **kwargs):
        super(YeastExperiment, self).__init__(*args, **kwargs)

    def _get_proteins_and_info(self):
        """Reads the ORFs and their sequences from the endpoint."""
        ans = self.query("""
        SELECT ?orf ?seq
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                 ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                 ocelot:sgd_id_has_sequence ?seq .
        }}
        """)
        assert ans and len(ans) and "results" in ans
        p_to_seq = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p = bindings[u"orf"].split(".")[-1]
            p_to_seq[p] = bindings[u"seq"]

        ans = self.query("""
        SELECT ?orf ?fun
        FROM <{default_graph}>
        WHERE {{
            ?orf a ocelot:sgd_id ;
                ocelot:sgd_id_has_type ocelot:sgd_feature_type.ORF ;
                ocelot:sgd_id_has_goslim ?fun .
        }}
        """)
        assert ans and len(ans) and "results" in ans
        p_to_fun = {}
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p = bindings[u"orf"].split(".")[-1]
            fun = bindings[u"fun"].split("#")[-1]
            if not p in p_to_fun:
                p_to_fun[p] = set([ fun ])
            else:
                p_to_fun[p].update([ fun ])

        p_to_feat = {}
        return set(p_to_seq.keys()), p_to_seq, p_to_fun, p_to_feat

    def _get_protein_interactions(self, manual_only = False):
        query = """
        SELECT ?orf1 ?orf2
        FROM <{default_graph}>
        WHERE {{
            ?orf1 a ocelot:sgd_id ;
                owl:sameAs ?feat1 .
            ?orf2 a ocelot:sgd_id ;
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
        print query
        ans = self.query(query)
        assert ans and len(ans) and "results" in ans
        pp_interactions = set()
        for bindings in ans["results"]["bindings"]:
            bindings = { k: self.cast(v) for k, v in bindings.items() }
            assert len(bindings) == 2
            p1 = bindings[u"orf1"].split(".")[-1]
            p2 = bindings[u"orf2"].split(".")[-1]
            pp_interactions.update([ (p1,p2), (p2,p1) ])
        return pp_interactions

    def run(self):
        """Run the yeast prediction experiment."""
        from pprint import pprint
        ps, p_to_seq, p_to_fun, p_to_feat = self._get_proteins_and_info()
        print _cls(self), ": found {} ps".format(len(ps))

        pp_interactions = self._get_protein_interactions(manual_only = True)
        density = float(len(pp_interactions)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} HQ pp interactions (positive density = {})".format(len(pp_interactions), density)

        pp_interactions = self._get_protein_interactions(manual_only = False)
        density = float(len(pp_interactions)) / (len(ps) * (len(ps) - 1))
        print _cls(self), ": found {} LQ pp interactions (positive density = {})".format(len(pp_interactions), density)

    # figure out the DDI-RRI (only use instance information)
    # match the PPI and DDI-RRI
    # compute composition features
    # **compute complexity features
    # **compute conservation (profile) features
    # **compute secondary features
    # compute subloc features (binary vec)
    # compute COG features (binary vec)
    # compute cell-cycle gene expression features (correlations)
    # compute environment-response gene expression features (correlations)
    # read in the Y2H raw data
    # read in the TAP-MS raw data
    # read in the Yip folds
    # write the SBR output (individuals, kernels, predicates, rules)
    # TODO: filter away short proteins
    # TODO: redundancy-reduce the protein dataset, map proteins to their
    #       representatives
    # TODO: form the CV folds
    # TODO: put all duplicates back into the training set
