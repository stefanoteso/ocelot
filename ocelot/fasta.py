# -*- coding: utf-8 -*-

from ocelot.utils import run_binary

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
