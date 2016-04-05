# -*- coding: utf-8 -*-

import os, errno, hashlib
import numpy as np
import multiprocessing as mp
import itertools as it
from ocelot.sequences import read_fasta, write_fasta, AMINOACIDS

# Taken from the string kernels package
BACKGROUND_AA_FREQ = (
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
)

def iterate_csv(path, num_skip = 0, **kwargs):
    """Thin wrapper around ``csv.DictReader`` with a couple more options."""
    import csv
    with open(path, "rU") as fp:
        for row_as_dict in csv.DictReader(fp, **kwargs):
            if num_skip > 0:
                num_skip -= 1
                continue
            yield row_as_dict

class InterProTSV(object):
    """Reader for `InterPro <http://www.ebi.ac.uk/interpro/>`_ tab-separated values files.

    The TSV files can be obtained with ``iprscan5_*.py`` or similar tools,
    which can be found `here <http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest>`_.

    .. todo::

        We should really integrate the ``iprscan5_*.py`` functionality here.
    """
    def read(self, path, allowed_sources = None):
        """
        :param path: path to the InterPro TSV file.
        :param allowed_sources: list of allowed domain sources (default: all).

        :returns: a list of ``(domain_id, evalue)`` pairs.

        Note that some domains may have ``evalue`` set to ``None``; this
        happens whenever InterPro does not provide a value.
        """
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

class Binary(object):
    """A simple wrapper around binary executables."""
    def __init__(self, path):
        self.path = path

    def run(self, args, shell = True):
        import subprocess
        pipe = subprocess.Popen(self.path + " " + " ".join(args),
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE,
                                shell = shell)
        out, err = pipe.communicate()
        ret = pipe.wait()
        return ret, out, err

# The amino acid alphabet (sorted as in PSSM files).
AMINOACIDS_PSSM = ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                   "M", "F", "P", "S", "T", "W", "Y", "V")

_PSSM_TO_LEX = {
    i: AMINOACIDS.index(s) for i, s in enumerate(AMINOACIDS_PSSM)
}

class PSSM(object):
    """A container for PSI-Blast PSSM profiles.

    :param targets: collection of targets (default: all)
    :param prior: list of priors on amino acid probabilities (default:
        BACKGROUND_AA_FREQ)
    :param gamma: amount of smoothing-by-prior (default: 0.8).
    """
    def __init__(self, **kwargs):
        TARGETS = ["residue", "score", "condp", "nlog_condp"]
        targets = kwargs.get("targets", TARGETS)
        if any(map(lambda target: not target in TARGETS, targets)):
            raise ValueError("all targets must in '{}'".format(TARGETS))
        self._targets = targets
        priors = kwargs.get("priors", BACKGROUND_AA_FREQ)
        if len(priors) != 20:
            raise ValueError("len(priors) must be 20")
        if any(map(lambda prior: not (0.0 <= prior <= 1.0), priors)):
            raise ValueError("priors must be in [0,1]")
        self._priors = priors
        gamma = kwargs.get("gamma", 0.8)
        if not (0.0 <= gamma <= 1.0):
            raise ValueError("gamma must be in [0,1]")
        self._gamma = gamma

    def read(self, path):
        """Reads a PSSM file.

        :param path: path to the file to be read.
        :returns: a list of lists, with one inner list for each PSSM residue.

        The output probabilities/scores are ordered as in ``AMINOACIDS``.
        """
        with open(path, "rt") as fp:
            lines = fp.readlines()[3:-6]
        infos = []
        old_position = None
        for line in lines:
            # Unfortunately, PSSM is one of those awful textual file formats
            # where the fields can not be extracted with a simple split()...
            position = int(line[0:5].strip()) - 1
            if old_position:
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
            info = []
            if "residue" in self._targets:
                info.append(residue)
            if "score" in self._targets:
                info.append(scores)
            if "condp" in self._targets:
                info.append(condps)
            if "nlog_condp" in self._targets:
                info.append(nlog_condps)
            infos.append(info)
        return infos

class PsiBlast(object):
    """A simple wrapper around NCBI PSI-Blast.

    The Blast database should be located in the path specified by either the
    ``BLASTDB`` environment variable or by the ``~/.ncbirc`` file (see
    `here <http://www.ncbi.nlm.nih.gov/books/NBK52640/>`_ for more
    information).

    :param cachedir: path of the cache directory (default: ``"/tmp/ocelot/psiblast"``)
    :param db: database name, e.g. `nr`.
    :param evalue: E-value as a string.
    :param matrix: similarity matrix to use.
    :param num_iterations: number of BLAST iterations, minimum 2.
    :param num_threads: number of CPU threads to use.

    The other ``kwargs`` are passed to the PSSM constructor.
    """
    def __init__(self, **kwargs):
        self._kwargs = kwargs

    def run(self, sequence):
        cachedir = self._kwargs.get("cachedir", "/tmp/ocelot/psiblast")
        try:
            os.makedirs(cachedir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise RuntimeError("can not create cache directory '{}': {}" \
                                    .format(cachedir, e))

        db              = self._kwargs.get("db", "nr")
        evalue          = self._kwargs.get("evalue", "10.0")
        matrix          = self._kwargs.get("matrix", "BLOSUM62")
        num_iterations  = self._kwargs.get("num_iterations", "2")
        num_threads     = self._kwargs.get("num_threads", str(mp.cpu_count()))
    
        m = hashlib.md5()
        m.update("_".join([sequence, db, evalue, matrix, num_iterations]))
        basename = os.path.join(cachedir, m.hexdigest())

        pssm_parser = PSSM()

        try:
            with open(basename + ".f", "rb") as fp:
                contents = map(str.strip, fp)
                if len(contents) != 1 or contents[0] != sequence:
                    raise ValueError("content mismatch")
            return pssm_parser.read(basename + ".ascii-pssm")
        except:
            pass

        with open(basename + ".f", "wb") as fp:
            fp.write(">temp\n{}\n".format(sequence))
        args = [
            "-query {}".format(basename + ".f"),
            "-out {}".format(basename + ".f.out"),
            "-out_pssm {}".format(basename + ".pssm"),
            "-out_ascii_pssm {}".format(basename + ".ascii-pssm"),
            "-db {}".format(db),
            "-evalue {}".format(evalue),
            "-matrix {}".format(matrix),
            "-num_iterations {}".format(num_iterations),
            "-num_threads {}".format(num_threads),
        ]
        ret, out, err = Binary("psiblast").run(args)
        if ret != 0:
            raise RuntimeError("psiblast exited with error '{}': {}".format(ret, err))
        return pssm_parser.read(basename + ".ascii-pssm")

class _TestPsiBlast(object):
    def test_output(self):
        psiblast = PsiBlast()
        psiblast.run("".join(AMINOACIDS))
