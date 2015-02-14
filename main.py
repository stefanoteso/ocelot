#!/usr/bin/env python2
# -*- coding: utf8 -*-

import os
import ocelot

from collections import namedtuple

def _checkdir(path):
    """Checks if a path identifies a directory."""
    return path and os.path.isdir(os.path.abspath(path))

def _filter_targets(ALL_TARGETS, ids):
    """Returns the targets identified by `ids`."""
    if ids == None:
        ids = ALL_TARGETS.keys()
    for id_ in ids:
        if not id_ in ALL_TARGETS:
            raise ValueError("invalid target id '{}'".format(id_))
    return [(id_, ALL_TARGETS[id_]) for id_ in ids]

def _make_rdf(args):
    """Builds the RDF dataset.

    :param args.src: path to the source data (i.e. root of all databases).
    :param args.dst: path to the directory to put the RDF files in.
    :param args.targets: list of targets to convert.
    :param args.force_update: force updating the RDF files even if they already
                              exist.
    """
    if not _checkdir(args.src):
        raise ValueError("Source directory '{}' does not exist".format(args.src))
    if not _checkdir(args.dst):
        raise ValueError("Target directory '{}' does not exist".format(args.dst))

    Target = namedtuple("Target", ("Converter", "kwargs"))
    ALL_TARGETS = {
        # DATABASES
        "sgd"       : Target(ocelot.SGDConverter, {}),
        "ipfam"     : Target(ocelot.IPfamConverter, {}),
        "sifts"     : Target(ocelot.SIFTSConverter, {}),
        "pdb"       : Target(ocelot.PDBConverter, {}),
        "biogrid"   : Target(ocelot.BioGRIDConverter, {
                        "path": "BioGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.2.112.mitab.txt"
                      }),
        "mint"      : Target(ocelot.PsiMiTabConverter, {
                        "path": "IMEx/IMEx-MINT-yeast.tab27"
                      }),
        "intact"    : Target(ocelot.PsiMiTabConverter, {
                        "path": "IMEx/IMEx-IntAct-yeast.tab27"
                      }),

        # DATASETS
        "yip09"     : Target(ocelot.Yip09Converter, {}),
    }
    for id_, target in _filter_targets(ALL_TARGETS, args.targets):
        target.kwargs["basename"] = id_
        converter = target.Converter(args.src, args.dst,
                                     force_update = args.force_update,
                                     **target.kwargs)
        converter.siphon()

def _upload_rdf(args):
    """Uploads the RDF data into a Virtuoso instance.

    It relies on the `isql` binary to upload the RDF files.

    :param args.src: path to the directory holding the RDF files.
    :param args.default_graph: default graph IRI.
    """
    from glob import glob

    if not _checkdir(args.src):
        raise ValueError("RDF directory '{}' does not exist".format(args.src))
    if args.default_graph == None:
        raise ValueError("upload-rdf requires a valid graph IRI")

    print "Uploading contents of '{}' into named graph <{}>".format(args.src, args.default_graph)

    turtles = []
    for path in glob(os.path.join(os.path.abspath(args.src), "*.ttl")):
        print "Adding '{}' to the upload list".format(path)
        turtles.append(path)

    isql_command = "log_enable(2);\n"
    for turtle in turtles:
        isql_command += "ttlp(file_open('{}'), '', '{}', 0);\n".format(
            turtle, args.default_graph)
    isql_command += "checkpoint;"

    with open("ocelot.isql", "wt") as fp:
        fp.write(isql_command)

    isql = ocelot.services.Binary("/opt/virtuoso/bin/isql")
    ret, out, err = isql.run("-U dba -P dba < ocelot.isql")
    if ret != 0:
        raise RuntimeError("upload-rdf failed, isql exited with '{}':\n{}\n{}\n".format(ret, out, err))

def _clear_rdf(args):
    """Clears a given graph URI.

    :param args.default_graph: default graph IRI.
    """
    if args.default_graph == None:
        raise ValueError("clear-rdf requires a valid graph IRI")

    with open("ocelot.isql", "wt") as fp:
        fp.write("SPARQL CLEAR GRAPH <{}>;".format(args.default_graph))

    isql = ocelot.services.Binary("/opt/virtuoso/bin/isql")
    ret, out, err = isql.run("-U dba -P dba < ocelot.isql")
    if ret != 0:
        raise RuntimeError("clear-rdf failed, isql exited with '{}':\n{}\n{}\n".format(ret, out, err))

def _run_experiment(args):
    """Converts the RDF data into a machine-learnable dataset.

    :param args.dst: destination directory for the results.
    :param args.subtargets: prediction targets.
    :param args.endpoint: Virtuoso SPARQL endpoint.
    :param args.default_graph: URI of the default graph.
    """
    Target = namedtuple("Target", ("Experiment", "kwargs"))
    ALL_TARGETS = {
        "yip"   : Target(ocelot.YipExperiment, {}),
        "yeast" : Target(ocelot.YeastExperiment, {}),
        "all"   : Target(ocelot.AllSpeciesExperiment, {}),
    }

    ALL_SUBTARGETS = {
        "ppi"   : None,
        "ddi"   : None,
        "rri"   : None,
        "pf-mf" : None,
        "pf-bp" : None,
        "pf-cc" : None,
    }
    subtargets = _filter_targets(ALL_TARGETS, args.subtargets)

    for id_, target in _filter_targets(ALL_TARGETS, args.targets):
        experiment = target.Experiment(args.src, args.dst,
                                       [st[0] for st in subtargets],
                                       endpoint = args.endpoint,
                                       default_graph = args.default_graph,
                                       **target.kwargs)
        experiment.run()

def main():
    import argparse as ap

    COMMANDS = {
        "make-rdf": _make_rdf,
        "upload-rdf": _upload_rdf,
        "clear-rdf": _clear_rdf,
        "run-experiment": _run_experiment,
    }

    parser = ap.ArgumentParser(description="Ocelot, the playful feline (also SBR for proteomics).")
    parser.add_argument("command", type=str, default=None,
                        help="What to do, any of {}".format(COMMANDS.keys()))
    parser.add_argument("-t", "--targets", nargs="+",
                        help="command-specific targets (default: all)")
    parser.add_argument("-y", "--subtargets", nargs="+",
                        help="target-specific sub-targets (default: all)")
    parser.add_argument("-s", "--src", type=str,
                        help="[make-rdf] where the input databases are")
    parser.add_argument("-d", "--dst", type=str,
                        help="[make-rdf] where to put the RDF files")
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="[make-rdf] force update the RDF files")
    parser.add_argument("-e", "--endpoint", type=str,
                        help="[upload-rdf,benchmark] set the SPARQL endpoint")
    parser.add_argument("-g", "--default-graph", type=str, default=None,
                        help="[upload-rdf] destination graph IRI")
    parser.add_argument("--debug", action="store_true",
                        help="[all] enable debugging output")
    args = parser.parse_args()

    command = COMMANDS.get(args.command)
    if not command:
        raise ValueError("invalid command '{}'".format)
    command(args)

if __name__ == "__main__":
    main()
