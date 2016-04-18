# -*- coding: utf-8 -*-

import pandas as pd
from ocelot.kernels import Kernel, SetKernel, SparseLinearKernel

def read_interpro(path, columns=None, allowed_sources=None):
    """Reads an InterPro [1]_ TSV file.

    TSV files can be obtained with the ``iprscan5_*.py`` scripts [2]_.

    Parameters
    ----------
    path : str
        Path to the InterPro TSV file.
    columns : list or None, optional. (defaults to None)
        Columns to keep. None equates to [DB, FAMILY, START, STOP, EVALUE,
        DATE, IPR_FAMILY, GO_TERMS, PATHWAYS].
    allowed_sources : list or None, optional. (defaults to None)
        List of allowed domain providers. None stands for all.

    Returns
    -------
    hits : pandas.DataFrame
        Each row is an InterPro hit. The EVALUE is set to None when the TSV
        file has no evalue.

    References
    ----------
    .. [1] `<http://www.ebi.ac.uk/interpro/>`_
    .. [2] `<http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest>`_
    """
    NAMES = [
        "QUERY_ID", "?2", "?3", "DB", "FAMILY", "DESCRIPTION", "START", "STOP",
        "EVALUE", "?10", "DATE", "IPR_FAMILY", "SHORT_DESCRIPTION", "GO_TERMS",
        "PATHWAYS"
    ]
    DEFAULT_COLUMNS = [
        "DB", "FAMILY", "START", "STOP", "EVALUE", "DATE", "IPR_FAMILY",
        "GO_TERMS", "PATHWAYS"
    ]

    df = pd.read_csv(path, sep="\t", names=NAMES)
    df = df[DEFAULT_COLUMNS if columns is None else columns]
    if allowed_sources is not None:
        df = df[df["SOURCE_DB"].isin(allowed_sources)]
    return df

class InterProKernel(Kernel):
    """A simple domain kernel built around InterPro.

    Parameters
    ----------
    ps : collection
        Ordered list of protein IDs.
    path : str
        Path to the directory holding the InterPro files.
    mode : str, optional
        One of "match", "count", "evalue", defaults to "match".
    allowed_sources : collection or None, optional
        List of allowed domain providers, defaults to all of them.
    default_score : float
        Score to use when no E-value is provided, defaults to 1.0.
    """
    def __init__(self, ps, path, mode="match", allowed_sources=None,
                 default_score=1.0, *args, **kwargs):
        if not mode in ("match", "count", "evalue"):
            raise ValueError("invalid mode '{}'".format(mode))
        self._path = path
        self._mode = mode
        self._allowed_sources = allowed_sources
        self._default_score = default_score
        super(InterProKernel, self).__init__(ps, *args, **kwargs)

    def _to_score(self, evalue):
        if evalue is None or evalue <= 0.0:
            return self._default_score
        return -np.log(evalue)

    def _compute_all(self):
        all_hits, num_missing = [], 0
        for p in self._entities:
            try:
                path = join(self._path, "{}.tsv.txt".format(p))
                df = read_interpro(path, self._allowed_sources)
            except IOError, e:
                domain_to_evalue = {}
                num_missing += 1
            raise NotImplementedError()

            if self._mode == "match":
                hits = set(domain_to_evalue.keys())
            elif self._mode == "count":
                hits = dict(Counter(domain_to_evalue.keys()))
            elif self._mode == "evalue":
                hits = {domain: self._to_score(evalue)
                        for domain, evalue in domain_to_evalue.iteritems()}
            all_hits.append(hits)

        if num_missing > 0:
            print "no interpro domains for '{}/{}' proteins" \
                    .format(num_missing, len(self))

        if self._mode == "match":
            return SetKernel(all_hits).compute()
        else:
            return SparseLinearKernel(all_hits).compute()
