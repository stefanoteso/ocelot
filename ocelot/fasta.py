# -*- coding: utf-8 -*-

def read_fasta(self, path):
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
            fp.write("{}\n{}\n".format(header.lstrip(">"), sequence))
