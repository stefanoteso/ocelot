# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.stats import entropy
from os.path import join
import multiprocessing as mp
import hashlib
from ocelot.utils import quietmkdir, validate, run_binary


# TODO add B (Asparagine or Aspartic Acid), Z (Glutamine or Glutamic Acid), J
# (Leucine or Isoleucine), U (Selenocystein), O (Pyrrolisine)
AA_INFO = pd.read_csv("data/aainfo.tsv", sep="	").set_index(u"aa")
"""Amino acid information, taken from [WikipediaAAa]_, [WikipediaAAb]_,
[Swart04]_."""

AA1_TO_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}
AA3_TO_AA1 = {aa3: aa1 for aa1, aa3 in AA1_TO_AA3.iteritems()}

AMINOACIDS = list(AA_INFO.index)
"""The amino acids alphabet, sorted lexicographically, one-letter format."""

AMINOACIDS3 = [AA1_TO_AA3[aa] for aa in AMINOACIDS]
"""The amino acids alphabet, sorted lexicographically, three-letter format."""

AMINOACIDS_PSSM = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                   "M", "F", "P", "S", "T", "W", "Y", "V"]
"""The amino acids alphabet, sorted as in PSSM files."""

_PSSM_TO_LEX = {
    i: AMINOACIDS.index(s) for i, s in enumerate(AMINOACIDS_PSSM)
}

AA_INFO.loc["?"] = AA_INFO.mean()

BACKGROUND_AA_FREQ = [
    0.0799912015849807, # A
    0.0171846021407367, # C
    0.0578891399707563, # D
    0.0638169929675978, # E
    0.0396348024787477, # F
    0.0760659374742852, # G
    0.0223465499452473, # H
    0.0550905793661343, # I
    0.060458245507428,  # K
    0.0866897071203864, # L
    0.0215379186368154, # M
    0.044293531582512,  # N
    0.0465746314476874, # P
    0.0380578923048682, # Q
    0.0484482507611578, # R
    0.0630028230885602, # S
    0.0580394726014824, # T
    0.0700241481678408, # V
    0.0144991866213453, # W
    0.03635438623143,   # Y
]

def read_fasta(path):
    """Reads a FASTA file (generator).

    Start-of-header symbols and end-of-sequence symbols (a.k.a. "*") are
    stripped.

    Parameters
    ----------
    path : str
        Path to the FASTA file.

    Returns
    -------
    data : list
        List of pairs of the form ``(header, sequence)``.
    """
    with open(path, "rt") as fp:
        header = None
        sequence = None
        for line in fp:
            line = line.strip()
            if len(line) == 0:
                continue
            elif line[0] == ">":
                if not header is None:
                    assert(sequence)
                    yield header, sequence.rstrip("*")
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        yield header, sequence.rstrip("*")

def write_fasta(path, data):
    """Writes a FASTA file.

    Parameters
    ----------
    path : str
        Path to the FASTA file.
    data : list
        List of pairs of the form ``(header, sequence)``.
    """
    with open(path, "wt") as fp:
        for header, sequence in data:
            if not header.startswith(">"):
                header = ">" + header
            fp.write("{}\n{}\n".format(header, sequence))

def get_composition_phi(sequence, what=None):
    """Composition features based on ``AA_INFO``.

    Parameters
    ----------
    sequence : str
        The amino acid sequences.
    what : list, optional. (defaults to ``None``)
        Target composition information. ``None`` means all of them.

    Returns
    -------
    phi : np.ndarray of shape (nres, nwhat)
        Array with one residue per row, one target feature per column.
    """
    if what is None:
        what = list(AA_INFO.columns)
    aas = [aa if aa in AA_INFO.index else "?" for aa in sequence.upper()]
    return AA_INFO.loc[aas][what].values

def get_low_complexity_phi(sequence, windows=None, epsilon=1e-10):
    """Low-complexity region features.

    Low-complexity regions are correlated with packing density and
    hydrophobicity and with the interface propensity of residues.
    Some useful references are [Wootton94]_, [Liao05]_, [Coletta10]_.

    Measured by computing the Shannon entropy over one or more windows of
    residues.

    Parameters
    ----------
    sequence :
        The amino acid sequence.
    windows : list, optional. (defaults to [8, ..., 128])
        WRITEME
    epsilon : float, optional. (defaults to 1e-10)
        WRITEME

    Returns
    -------
    phi : np.ndarray of shape (nres, nwindows)
        Array with one residue per row, one window per column.
    """
    if windows is None:
        windows = [8, 16, 32, 64, 128]
    phi = []
    for window in windows:
        window_phi = []

        padding = "?" * window
        s = padding + sequence + padding
        for i in range(len(sequence)):
            i_min = len(padding) + i - window
            i_max = len(padding) + i + window + 1
            w = s[i_min:i_max]

            counts = [len(filter(lambda res: res == aa, w)) + epsilon
                      for aa in AMINOACIDS]

            h = entropy(np.array(counts, dtype=np.float32) / sum(counts))
            if np.abs(h) < epsilon:
                h = 0.0
            window_phi.append(h)

        phi.append(window_phi)
    return np.array(phi).T

def get_sequence_phi(get_phi, sequences, *args, **kwargs):
    """Per-sequence average features.

    Parameters
    ----------
    get_phi : callable
        The feature construction function.
    sequences : list
        List of amino acid sequences.
    All other parameters are passed to get_phi().

    Returns
    -------
    phi : np.ndarray of shape (nsequences, nfeatures)
        The per-sequence features, computed by averaging the per-residue
        features.
    """
    return np.array([np.mean(get_phi(sequence, *args, **kwargs), axis=0)
                     for sequence in sequences])

def _read_cdhit_results(path):
    with open(path + ".clstr") as fp:
        clusters, members = [], None
        for line in fp:
            if line.startswith(">"):
                if not members is None:
                    clusters.append(members)
                members = set()
            else:
                words = line.split()
                if len(words) == 4:
                    _, _, header, _ = words
                    perc = 1.0
                elif len(words) == 5:
                    _, _, header, _, perc = words
                    perc = float(perc.strip("%")) / 100.0
                else:
                    raise SyntaxError("unexpected line '{}'".format(line))
                header = header[1:-3]
                members.add((header, perc))
        clusters.append(members)
    return [header for header, _ in read_fasta("temp.cdhit")], clusters

def cdhit(pairs, threshold=0.8, cdhit_path="/usr/bin/cdhit"):
    """Cluster protein sequences using CD-HIT [1]_.

    Parameters
    ----------
    pairs : list
        List of pairs of the form (id, sequence)
    threshold : float, optional. (defaults to 0.8)
        CD-HIT clustering threshold.
    " ".cdhit_path : str, optional. (defaults to ``/usr/bin/cdhit``)
        Path to the cdhit binary.

    Returns
    -------
    exemplars : list
        WRITEME
    clusters : list
        WRITEME

    References
    ----------
    .. [1] `cdhit <http://weizhongli-lab.org/cd-hit/>`_
    """
    # TODO use tempfile, get rid of it
    write_fasta("temp.fasta", pairs)

    # TODO add support for the other options
    args = [ "-i {}".format("temp.fasta"),
             "-o {}".format("temp.cdhit"),
             "-c {}".format(threshold) ]

    ret, out, err = run_binary(cdhit_path, args)
    if ret != 0:
        raise RuntimeError("cdhit exited with errno '{}'".format(ret))

    return _read_cdhit_results("temp.cdhit")

def read_pssm(path, columns=None, priors=None):
    """Reads a PSSM (position-specific scoring matrix) file.

    Parameters
    ----------
    path : str
        Path to the PSSM file.
    columns : list or None. (defaults to None)
        List of columns to keep.
    priors : numpy.ndarray of shape (20,) or None. (defaults to None)
        Amino acid background frequencies. None means BACKGROUND_AA_FREQ.
    gamma : float. (defaults to 0.8)
        Amount of smoothing.

    Returns
    -------
    pssm : pandas.DataFrame
        A DataFrame representing the PSSM data.
    """
    NAMES = ["residue", "score", "condprob", "nlog_condprob"]

    columns = validate(NAMES, columns)
    if priors is None:
        priors = BACKGROUND_AA_FREQ
    if len(priors) != 20:
        raise ValueError("len(priors) must be 20")
    if not all(0.0 <= prior <= 1.0 for prior in priors):
        raise ValueError("priors must be in [0,1]")
    if not 0.0 <= gamma <= 1.0:
        raise ValueError("gamma must be in [0,1]")

    with open(path, "rt") as fp:
        lines = fp.readlines()[3:-6]

    infos, old_position = [], None
    for line in lines:
        position = int(line[0:5].strip()) - 1
        if old_position is not None:
            assert position == (old_position + 1), "watch your PSSMs"
        old_position = position

        residue = line[5:7].strip()

        scores, condps, nlog_condps = \
            [ 0. ] * 20, [ 0. ] * 20, [ 0. ] * 20
        for i in xrange(20):
            score = float(line[9+i*3:9+i*3+3].strip())
            condp = float(line[70+i*4:70+i*4+4].strip()) / 100.0
            condp = self._gamma * condp + (1-self._gamma) * self._priors[i]
            assert 0.0 <= condp <= 1.0
            if condp == 0.0:
                print "Warning: zero conditional probability found in PSSM"
                nlog_condp = 1e13 # supposedly a **huge** value
            else:
                nlog_condp = -np.log(condp)
            assert nlog_condp >= 0.0
            scores[_PSSM_TO_LEX[i]] = score
            condps[_PSSM_TO_LEX[i]] = condp
            nlog_condps[_PSSM_TO_LEX[i]] = nlog_condp

        infos.append([residue, scores, condps, nlog_condps])

    return pd.DataFrame(infos, columns=NAMES)[columns]

def run_psiblast(sequence, db="nr", evalue=10.0, matrix="BLOSUM62",
                 num_iters=2, num_threads=0, cache="/tmp/psiblast",
                 psiblast="/usr/bin/psiblast"):
    """Runs NCBI PSI-Blast.

    Parameters
    ----------
    sequence : str
        Sequence to align.
    db : str, optional. (defaults to "nr")
        Name of or path to the BLAST database.
    evalue : float, optional. (defaults to 10)
        E-value.
    matrix : str, optional. (defaults to "BLOSUM62")
        Similarity matrix.
    num_iters : int, optional. (defaults to 2)
        Number of BLAST iterations, minimum 2.
    num_threads : int, optional. (defaults to 0)
        Number of threads to use. Non-positive means all.
    cache : str, optional. (defaults to "/tmp/psiblast")
        Path to temporary directory.
    psiblast : str, optional. (defaults to "/usr/bin/psiblast")
        Path to the psiblast binary.

    Returns
    -------
    pssm :
        WRITEME
    """
    quietmkdir(cache)

    h = hashlib.md5()
    h.update("_".join([sequence, db, evalue, matrix, num_iters]))
    basename = join(cache, h.hexdigest())

    write_fasta(basename + ".f", [("temp", sequence)])

    if num_threads <= 0:
        num_threads = mp.count_cpus()

    args = [
        "-query {}".format(basename + ".f"),
        "-out {}".format(basename + ".f.out"),
        "-out_pssm {}".format(basename + ".pssm"),
        "-out_ascii_pssm {}".format(basename + ".ascii-pssm"),
        "-db {}".format(db),
        "-evalue {}".format(evalue),
        "-matrix {}".format(matrix),
        "-num_iterations {}".format(num_iters),
        "-num_threads {}".format(num_threads),
    ]

    ret, out, err = Binary("psiblast").run(args)
    if ret != 0:
        raise RuntimeError("psiblast exited with error '{}': {}".format(ret, err))

    return read_pssm(basename + ".ascii-pssm")

def run_blastall(sequences, mkblast="/usr/bin/mkblast", blastp="/usr/bin/blastp"):

    # ncbi-blast-2.2.24+/bin/makeblastdb -in good_proteins.fasta -dbtype prot -out my_prot_blast_db
    args = [
        "-in sequences.fasta",
        "-dbtype prot",
        "-out sequences",
    ]
    ret, out, err = run_binary(mkblast, args)
    if ret != 0:
        raise RuntimeError("mkblast exited with errno '{}'".format(ret))

    args = [
        "-db {}".format(db),
        "-query {}".format("sequences.fasta"),
        "-outfmt {}".format(6),
        "-out".format("allvsall.tsv"),
        "-num_threads {}".format(num_threads),
    ]
    ret, out, err = run_binary(blastp, args)
    if ret != 0:
        raise RuntimeError("mkblast exited with errno '{}'".format(ret))
