# -*- coding: utf-8 -*-

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
    """A simple wrapper around fasta files."""
    def read(self, path):
        """Reads the contents of a fasta file.

        It generates pairs of the form `header, sequence`.

        It strips `*` end-of-sequence symbols from the end of the
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
        """Writes a fasta file out of a `header->sequence` dictionary."""
        with open(path, "wt") as fp:
            for header, sequence in data:
                fp.write("{}\n{}\n".format(header, sequence))

class PSSM(object):
    """A container for PSI-Blast PSSM profiles."""
    def __init__(self, path):
        pos_to_info = {}
        with open(path, "rt") as fp:
            for line in fp:
                print line.rstrip()
                words = line.rstrip().split()
                if len(words) != 44:
                    continue
                position = int(words[0]) - 1
                pos_to_info[position] = {
                    "res"   : words[1],
                    "data"  : words[2:2+20],
                    "info"  : words[2+20:2+40],
                    "unk1"  : float(words[2+40]),
                    "unk2"  : float(words[2+41]),
                }
        self.pos_to_info = pos_to_info
    def __getitem__(self, key):
        if isinstance(key, int):
            self._pos_to_info[key]
        elif isinstance(key, slice):
            for i in xrange(key.start, key.stop, key.step):
                self._pos_to_info[key]
        else:
            raise TypeError, "invalid key type '{}'".format(type(key))

class PCL(object):
    """Reads a Stanford PCL gene expression file.

    The file is simply a TSV where the first three columns are fixed,
    followed by a variable number of columns (one per condition).

    :param path: path to the PCL file.

    XXX we use the NAME column rather than the YORF column, as the NAME's
    are unique in our files while the YORFs are not. No idea why, really.

    XXX we also ignore the GWEIGHT, its value seems rather arbitrary
    anyway.

    *References*

    .. [PCL] http://smd.princeton.edu/help/formats.shtml#pcl
    """
    def read(self, path):
        import numpy as np

        FIELDS = ("ORF", "NAME", "GWEIGHT")
        orf_to_expression = {}
        num_conditions = -1
        ln = -1
        for row in iterate_csv(path, delimiter = "\t",
                               fieldnames = FIELDS):
            ln += 1
            if ln < 2:
                continue

            orf = row["NAME"]
            assert not orf in orf_to_expression, orf

            expression_levels = map(float, row[None])
            if num_conditions < 0:
                num_conditions = len(expression_levels)
            else:
                assert num_conditions == len(expression_levels)

            orf_to_expression[orf] = np.array(expression_levels)
        return orf_to_expression, num_conditions

class PsiBlast(object):
    """A simple wrapper around NCBI PSI-Blast.

    NOTE: the database should be located in the path specified by either the
    `BLASTDB` environment variable or by the `~/.ncbirc` INI file[1].

    :param db: database name, e.g. `nr`.
    :param evalue: E-value as a string.
    :param matrix: similarity matrix to use.
    :param num_iterations: number of BLAST iterations, minimum 2.
    :param num_threads: number of CPU threads to use.

    *References*

    [1] http://www.ncbi.nlm.nih.gov/books/NBK52640/
    """
    def __init__(self, cache = "/tmp/psiblast", **kwargs):
        for k in kwargs:
            assert k in ("db", "evalue", "matrix", "num_iterations", "num_threads")
        self.cache = cache
        self.kwargs = kwargs
    def _compute_one(self, sequence):
        """Computes the PSSM profile of a protein sequence.

        TODO: memoize (or cache).

        :param sequence: the protein sequence.
        :returns: a PSSM object.
        """
        import os, multiprocessing as mp, hashlib

        db              = self.kwargs.get("db", "nr")
        evalue          = self.kwargs.get("evalue", "10.0")
        matrix          = self.kwargs.get("matrix", "BLOSUM62")
        num_iterations  = self.kwargs.get("num_iterations", "2")
        num_threads     = self.kwargs.get("num_threads", str(mp.cpu_count()))

        m = hashlib.md5()
        m.update("_".join([db, evalue, matrix, num_iterations,
                           sequence]))
        basename = os.path.join(self.cache, m.hexdigest())
        try:
            with open(basename + ".in", "r") as fp:
                stored_sequence = map(str.strip, fp.readlines())[0]
                if sequence != stored_sequence.strip():
                    raise RuntimeError, "content mismatch"
            return PSSM(basename + ".pssm")
        except:
            pass
        try:
            with open(basename + ".in", "w") as fp:
                fp.write(">temp\n")
                fp.write(sequence)
            args = [ "-query {}".format(basename + ".in"),
                     "-out {}".format(basename + ".out"),
                     "-out_ascii_pssm {}".format(basename + ".pssm"),
                     "-db {}".format(db),
                     "-evalue {}".format(evalue),
                     "-matrix {}".format(matrix),
                     "-num_iterations {}".format(num_iterations),
                     "-num_threads {}".format(num_threads) ]
            ret, out, err = Binary("psiblast").run(args)
            assert ret == 0
        except Exception, e:
            print e
        return PSSM(basename + ".pssm")

class BlastSeqSim(object):
    """
    Uses blastall to compute pairwise sequence similarity.
    """

    def __init__(self, db, binpath = "/usr/bin/legacy_blast blastall", **kwargs):
        self.blastall = Binary(binpath)
        self.db = db

    def run(self, sequences):
        fasta = FASTA()
        fasta.write("temp.fasta", sequences)

#class DSSP(object):
#    """Handles DSSP services.
#
#    Provides a facility for computing DSSP files out of PDB and mmCIF
#    files, whether gzip'ed or not, and for parsing the resulting DSSP
#    files.
#
#    For more information on the various fields, see [1].
#
#    Tested with DSSP 2.2.1.
#
#    [1] http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html
#    """
#
#    def __init__(self, binpath = "/usr/bin/dssp"):
#        self.binary = Binary(binpath)
#
#    def parse(self, fp):
#        @unique
#        class C(Enum):
#            # DSSP internal line/record number
#            NUM    = (0, 5,    lambda s: int(s.strip()))
#            # DSSP residue number (use this); note that it may not
#            # start from either 0 nor 1...
#            RESNUM    = (5, 10,    lambda s: int(s.strip()))
#            # PDB insertion code
#            INSCODE    = (10, 11,    lambda s: s.strip())
#            # Protein chain
#            CHAIN    = (11, 12,    lambda s: s.strip())
#            # Residue amino-acid
#            AA    = (12, 14,    lambda s: s.strip())
#            # ?
#            STRUCT    = (14, 25,    lambda s: None)
#            # Residue number of the first bridge partner
#            BP1    = (25, 29,    lambda s: float(s.strip()))
#            # Residue number of the second bridge partner
#            BP2    = (29, 33,    lambda s: float(s.strip()))
#            # Sheet label
#            SHEET    = (34, 34,    lambda s: s.strip())
#            # Number of water molecules in contact with the residue
#            # _or_ residue water-exposed surface in square Angstrom
#            ACC    = (34, 38,    lambda s: float(s.strip()))
#            # ?
#            H_NHO1    = (38, 50,    lambda s: None)
#            H_OHN1    = (50, 61,    lambda s: None)
#            H_NHO2    = (61, 72,    lambda s: None)
#            H_OHN2    = (72, 83,    lambda s: None)
#            # Cosine of the angle between C=O of residue i and C=O
#            # of residue i-1. It is close to +1 for alpha helices
#            # and close to -1 for beta sheets.
#            TCO    = (83, 91,    lambda s: float(s.strip()))
#            # Virtual bond angle (bend angle) defined by the three
#            # carbon-alpha atoms of residues i-2, i, i+2
#            KAPPA    = (91, 97,    lambda s: float(s.strip()))
#            # Virtual torsion angle (dihedral angle) defined by the
#            # four carbon-alpha atoms of residues i-1, i, i+1, i+2
#            ALPHA    = (97, 103,    lambda s: float(s.strip()))
#            # IUPAC peptide backbone torsion angle
#            PHI    = (103, 109,    lambda s: float(s.strip()))
#            # IUPAC peptide backbone torsion angle
#            PSI    = (109, 115,    lambda s: float(s.strip()))
#            # Carbon-alpha coordinates
#            XCA    = (115, 122,    lambda s: float(s.strip()))
#            YCA    = (122, 129,    lambda s: float(s.strip()))
#            ZCA    = (129, 136,    lambda s: float(s.strip()))
#
#        data = {}
#
#        start = False
#        for line in fp:
#            if "#" in line:
#                start = True
#                continue
#            if not start:
#                continue
#            if "!*" in line:
#                # An end-of-chain marker
#                continue
#            line = line
#
#            row = {}
#            for c, (pos0, pos1, func) in [ (c, c.value) for c in C ]:
#                row[c] = func(line[pos0:pos1])
#
#            key = (row[C.CHAIN], row[C.RESNUM])
#            if not key in data:
#                data[key] = row
#
#        assert(start)
#
#        return data
#
#    def compute(self, in_path, out_path = None):
#
#        if not out_path:
#            out_path = in_path + ".dssp"
#
#        ret, out, err = self.binary.run("-i {} -o {}".format(in_path, out_path))
#        if ret != 0:
#            raise RuntimeError("DSSP: dssp exited with error code {}:\n\nstdout:\n{}\n\nstderr:\n{}\n\n".format(ret, out, err))
#
#        with open(out_path, "rb") as fp:
#            data = self.parse(fp)
#
#        return data

# class SignalP(object):
#     """Support for psort 3.1/4.0.
#     
#     Original author: Karunakar Bayyapu (karun AT cbs DOT dtu DOT dk), taken
#     from [0], version version 3.1/4.0 ws0, 2012-03-23.
#     
#     [0] http://www.cbs.dtu.dk/ws/SignalP4/examples/signalp_python
#     
#     [1] Petersen et al. "SignalP 4.0 - Discrimination between Signal
#         Peptides and Transmembrane Regions", 2011.
#     """
#     def __init__(self):
#         pass
#     def __getitem__(self, query):
#         pass
#     WSDL = 'http://www.cbs.dtu.dk/ws/SignalP4/SignalP4_4_0_ws0.wsdl'
# 
#     @unique
#     class ORGANISM(Enum):
#         EUKARIOTA    = "euk"
#         GRAM_NEGATIVE    = "gram-"
#         GRAM_POSITIVE    = "gram+"
# 
#     @unique
#     class METHOD(Enum):
#         HMM        = "best"
#         NEURALNET    = "notm"
# 
#     def __init__(self, **kwargs):
#         self.organism = kwargs.get("organism", self.ORGANISM.EUKARIOTA)
#         if not self.organism in self.ORGANISM:
#             raise ValueError("invalid organism '{}'".format(self.organism))
# 
#         self.method = kwargs.get("method", self.METHOD.HMM)
#         if not self.method in self.METHOD:
#             raise ValueError("invalid method '{}'".format(self.method))
# 
#         self.client = suds.client.Client(self.WSDL, cache = None)
#         self.id = 0
# 
#     def query(self, id, sequence):
#         print "DEBUG: querying '{}' as '{}'".format(sequence, id)
# 
#         client = self.client
# 
#         sequence = client.factory.create("runService.parameters.sequencedata.sequence")
#         sequence.id = "{}{}".format(id, self.id)
#         sequence.seq = sequence
#         self.id += 1
# 
#         pprint(sequence)
# 
#         request = client.factory.create("runService.parameters")
#         request.organism = self.organism.value
#         request.method = self.method.value
#         request.sequencedata.sequence = [ sequence ]
# 
#         pprint(request)
# 
#         client.service.runService(request)
# 
#         sys.exit(1)
# 
#         sequence = client.factory.create("pollQueue.job")
#         sequence.jobid = "whatever"
# 
#         response = self.client.service.pollQueue(sequence)
# 
# 
#         sequence = client.factory.create("fetchResult.job")
#         sequence.jobid = "whatever"
# 
#         response = self.client.service.fetchResult(sequence)
# 
# 
# 
# class ReProf(object):
#     pass
# 
#     def compute_ss_reprof(self, sequence):
# 
#         @unique
#         class C(IntEnum):
#             NO    = 0    # position, starting from 1
#             AA    = 1    # residue
# 
#             PHEL    = 2    # secondary structure (H = Helix, E = Sheet, L = Loop)
#             RI_S    = 3    # Reliability index, in [0,9] (low to high)
#             PH    = 4    # Probability of H, in [0,1]
#             PE    = 5    # Probability of E, in [0,1]
#             PL    = 6    # Probability of L, in [0,1]
# 
#             PACC    = 7    # Absolute solvent accessibility
#             PREL    = 8    # Relative solvent accessibility
#             P10    = 9    # Relative accessibility, in [0-9] (buried to exposed)
#             RI_A    = 10    # Reliability index, in [0,9] (low to high)
#             PBE    = 11    # Two-state accessibility (b = buried, e = exposed)
#             PBIE    = 12    # Three-state accessibility (b = buried, i = intermediate, e = exposed)
# 
#         path = self.os.path.abspath("temp.fasta")
#         with open("temp.fasta", "wt") as fp:
#             fp.write(">temp\n{}".format(sequence))
#             fp.flush()
# 
#         ret, out, err = self.call("reprof -i {} -o {}.reprof".format(path, path))
#         assert(ret == 0 and err == None)
# 
#         csv = ocelot.utils.CSV(comment = "#", columns = C)
# 
#         phis = []
#         for row in csv.read("temp.fasta.reprof"):
#             if row[C.NO] == "No":
#                 continue
#             phis.append([
#                 float(row[C.PH]) / 100.0,
#                 float(row[C.PE]) / 100.0,
#                 float(row[C.PL]) / 100.0,
#                 float(row[C.RI_S]) / 10.0,
#                 float(row[C.PACC]),
#                 float(row[C.PREL]),
#                 float(row[C.P10]) / 10.0,
#                 float(row[C.RI_A]) / 10.0,
#             ])
# 
#         return np.array(phis)
# 
