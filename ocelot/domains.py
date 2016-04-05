# -*- coding: utf-8 -*-

import pandas as pd

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
