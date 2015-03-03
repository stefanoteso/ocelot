# -*- coding: utf-8 -*-

import numpy as np

# The amino acid alphabet (sorted lexicographically).
AMINOACIDS = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
              "Q", "R", "S", "T", "V", "W", "Y")

# The amino acid alphabet (sorted as in PSSM files).
AMINOACIDS_PSSM = ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                   "M", "F", "P", "S", "T", "W", "Y", "V")

_PSSM_TO_LEX = {
    i: AMINOACIDS.index(s) for i, s in enumerate(AMINOACIDS_PSSM)
}

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

def _cls(obj):
    return type(obj).__name__

def iterate_csv(path, num_skip = 0, **kwargs):
    import csv
    with open(path, "rU") as fp:
        reader = csv.DictReader(fp, **kwargs)
        for row_as_dict in reader:
            if num_skip > 0:
                num_skip -= 1
                continue
            yield row_as_dict

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

class FASTA(object):
    """A simple wrapper around FASTA files."""
    def read(self, path):
        """Reads the contents of a FASTA file.

        :param path: path to the FASTA file to be read.
        :returns: pairs of the form ``(header, sequence)``, as a generator.

        Note that it strips ``*`` end-of-sequence symbols from the end of the
        sequence."""
        with open(path, "rt") as fp:
            header = None
            sequence = None
            for line in fp:
                line = line.strip()
                if len(line) == 0:
                    continue
                elif line[0] == ">":
                    if header != None:
                        assert(sequence)
                        yield header, sequence.rstrip("*")
                    header = line[1:]
                    sequence = ""
                else:
                    sequence += line
            yield header, sequence.rstrip("*")

    def write(self, path, data):
        """Writes a FASTA file.

        :param path: path to the FASTA file to be written.
        :param data: a list of ``(header, sequence)`` pairs.
        """
        with open(path, "wt") as fp:
            for header, sequence in data:
                fp.write(">{}\n{}\n".format(header.lstrip(">"), sequence))

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

class PCL(object):
    """Thin wrapper around `PCL <http://smd.princeton.edu/help/formats.shtml#pcl>`_
    gene expression files.

    The file is simply a TSV where the first three columns are fixed,
    followed by a variable number of columns (one per condition).
    """
    def read(self, path):
        """Reads a PCL file into a dictionary.

        .. note::

            This class uses the ``NAME`` column, rather than the ``YORF`` colum,
            as the key to the dictionary. The ``YORF``'s tend not to be unique
            within the file.

        .. note::

            This class currently ignores the ``GWEIGHT``.

        :param path: path to the PCL file.
        :returns: a {``name``: ``expression-levels``} dictionary.
        """
        FIELDS = ("ORF", "NAME", "GWEIGHT")
        orf_to_expression = {}
        num_conditions = -1
        for row in iterate_csv(path, delimiter = "\t", fieldnames = FIELDS,
                               num_skip = 2):
            orf = row["NAME"]
            assert not orf in orf_to_expression, orf

            expression_levels = map(float, row[None])
            if num_conditions < 0:
                num_conditions = len(expression_levels)
            else:
                assert num_conditions == len(expression_levels)

            orf_to_expression[orf] = np.array(expression_levels)
        return orf_to_expression, num_conditions

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

class CDHit(object):
    """A wrapper around the local `cdhit <http://weizhongli-lab.org/cd-hit/>`_
    installation.

    :param path: path to the ``cdhit`` binary (default: ``/usr/bin/cdhit``)

    .. todo:
        Use tempfile.

    .. todo:
        Get rid of the temporary files.

    .. todo:
        Add support for the gazillion missing options.
    """
    def __init__(self, path = "/usr/bin/cdhit"):
        self.cdhit = Binary(path)
    def run(self, pairs, **kwargs):
        """Cluster protein sequences using CD-HIT.

        :param pairs: pairs of the form ``(id, sequence)``.
        :param threshold: clustering threshold (default: ``0.9``).
        :returns: list of cluster representatives, list of clusters 
        """
        fasta = FASTA()
        fasta.write("temp.fasta", pairs)
        args = [ "-i {}".format("temp.fasta"),
                 "-o {}".format("temp.cdhit"),
                 "-c {}".format(kwargs.get("threshold", 0.8)) ]
        ret, out, err = self.cdhit.run(args)
        if ret != 0:
            raise RuntimeError("cdhit exited with errno '{}'".format(ret))

        with open("temp.cdhit.clstr") as fp:
            clusters = []
            members = None
            for line in fp:
                if line.startswith(">"):
                    if members != None:
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

        return [ k for k, v in fasta.read("temp.cdhit") ], clusters
