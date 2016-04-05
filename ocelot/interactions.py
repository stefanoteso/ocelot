# -*- coding: utf8 -*-

import pandas as pd

# NOTE all columns are mandatory.

_NAMES_25 = [
    "ID_INTERACTOR_A",
    "ID_INTERACTOR_B",
    "ALT_IDS_INTERACTOR_A",
    "ALT_IDS_INTERACTOR_B",
    "ALIASES_INTERACTOR_A",
    "ALIASES_INTERACTOR_B",
    "INTERACTION_DETECTION_METHOD",
    "PUBLICATION_1ST_AUTHOR",
    "PUBLICATION_IDENTIFIERS",
    "TAXID_INTERACTOR_A",
    "TAXID_INTERACTOR_B",
    "INTERACTION_TYPES",
    "SOURCE_DATABASE",
    "INTERACTION_IDENTIFIERS",
    "CONFIDENCE_VALUES",
]

_NAMES_26 = _NAMES_25 + [
    "COMPLEX_EXPANSION",
    "BIOLOGICAL_ROLE_A",
    "BIOLOGICAL_ROLE_B",
    "EXPERIMENTAL_ROLE_A",
    "EXPERIMENTAL_ROLE_B",
    "INTERACTOR_TYPE_A",
    "INTERACTOR_TYPE_B",
    "XREF_A",
    "XREF_B",
    "XREF_INTERACTION",
    "ANNOTATION_A",
    "ANNOTATION_B",
    "ANNOTATION_INTERACTION",
    "TAXID_HOST",
    "INTERACTION_PARAMETERS",
    "CURATION_CREATION_TIME",
    "CURATION_UPDATE_TIME",
    "CHECKSUM_A",
    "CHECKSUM_B",
    "CHECKSUM_INTERACTION",
    "IS_NEGATIVE",
]

_NAMES_27 = _NAMES_26 + [
    "FEATURES_A",
    "FEATURES_B",
    "STOICHIOMETRY_A",
    "STOICHIOMETRY_B",
    "PARTICIPANT_IDENTIFICATION_METHOD_A",
    "PARTICIPANT_IDENTIFICATION_METHOD_B",
]

def read_psimi(path, version=None):
    """Reads a PSI-MI TAB file.

    It supports the following versions of the PSI-MI TAB standard:
    - `2.5 <https://code.google.com/p/psicquic/wiki/MITAB25Format>`_
    - `2.6 <https://code.google.com/p/psicquic/wiki/MITAB26Format>`_
    - `2.7 <https://code.google.com/p/psicquic/wiki/MITAB27Format>`_

    Tested with PSI-MI TAB data obtained from BioGRID, IMEx/IntAct, IMEx/MINT,
    MIPS.

    Parameters
    ----------
    path : str
        Path to the PSI-MI TAB file.
    version : str or None, optional. (defaults to None)
        Version of the PSI-MI format to use. Can be "2.5", "2.6", or "2.7".
        None

    Returns
    -------
    data : pandas.DataFrame of shape (nproteins, values)
        The interaction data.
    version : str
        Actual version of the PSI-MI data.
    """
    df = pd.read_csv(path, sep="\t", header=0)
    if len(df.columns) == len(_NAMES_25):
        assert version is None or version == "2.5"
        df.columns = _NAMES_25
        version = "2.5"
    elif len(df.columns) == len(_NAMES_26):
        assert version is None or version == "2.6"
        df.columns = _NAMES_26
        version = "2.6"
    elif len(df.columns) == len(_NAMES_27):
        assert version is None or version == "2.7"
        df.columns = _NAMES_27
        version = "2.7"
    else:
        raise ValueError("invalid PSI-MI file")
    return df, version
