#!/usr/bin/env python2
# -*- coding: utf8 -*-

from os.path import join, abspath
import numpy as np
from collections import namedtuple
import ocelot
from ocelot.utils import run_binary
import experiments

def _make_rdf(args):
    """Builds the RDF dataset.

    :param args.src: path to the source data (i.e. root of all databases).
    :param args.dst: path to the directory to put the RDF files in.
    :param args.targets: list of targets to convert.
    :param args.force_update: force updating the RDF files even if they already
                              exist.
    """
    if not ocelot.checkdir(args.src):
        raise ValueError("Source directory '{}' does not exist".format(args.src))
    if not ocelot.checkdir(args.dst):
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
        "string"    : Target(ocelot.STRINGConverter, {
                        "taxon": "4932", "version": "9.1"
                      }),

        # DATASETS
        "yip09"     : Target(ocelot.Yip09Converter, {}),
        "cafa13"    : Target(ocelot.CAFA13Converter, {}),
    }
    for id_ in ocelot.validate(ALL_TARGETS.keys(), args.targets):
        target = ALL_TARGETS[id_]
        target.kwargs["basename"] = id_
        converter = target.Converter(args.src, args.dst,
                                     force_update = args.force_update,
                                     **target.kwargs)
        converter.siphon()

def _start_virtuoso(args):
    """Starts a Virtuoso instance"""
    section = ocelot.config["virtuoso"]
    virtuoso, ini = section[u"virtuoso"], section[u"ini"]
    print "About to start virtuoso (with '{}')".format(ini)
    ret, out, err = run_binary(virtuoso,["+configfile {}".format(ini)])
    if ret != 0:
        raise RuntimeError("virtuoso-t exited with error code '{}'".format(ret))

def _upload_rdf(args):
    """Uploads the RDF data into a Virtuoso instance.

    It relies on the `isql` binary to upload the RDF files.

    :param args.src: path to the directory holding the RDF files.
    :param args.default_graph: default graph IRI.
    """
    from glob import glob

    if not ocelot.checkdir(args.src):
        raise ValueError("RDF directory '{}' does not exist".format(args.src))
    if args.default_graph is None:
        raise ValueError("upload-rdf requires a valid graph IRI")

    print "Uploading contents of '{}' into named graph <{}>".format(args.src, args.default_graph)

    turtles = []
    for path in glob(join(abspath(args.src), "*.ttl")):
        print "Adding '{}' to the upload list".format(path)
        turtles.append(path)

    isql_command = "log_enable(2);\n"
    for turtle in turtles:
        isql_command += "ttlp(file_open('{}'), '', '{}', 0);\n".format(
            turtle, args.default_graph)
    isql_command += "checkpoint;"

    with open("ocelot.isql", "wt") as fp:
        fp.write(isql_command)

    ret, out, err = run_binary("isql", [ "-U dba -P dba < ocelot.isql" ])
    if ret != 0:
        raise RuntimeError("upload-rdf failed, isql exited with '{}':\n{}\n{}\n".format(ret, out, err))

def _clear_graph(args):
    """Clears a given graph URI.

    :param args.default_graph: default graph IRI.
    """
    if args.default_graph is None:
        raise ValueError("clear-rdf requires a valid graph IRI")

    with open("ocelot.isql", "wt") as fp:
        fp.write("SPARQL CLEAR GRAPH <{}>;".format(args.default_graph))

    ret, out, err = run_binary("isql", [ "-U dba -P dba < ocelot.isql" ])
    if ret != 0:
        raise RuntimeError("clear-rdf failed, isql exited with '{}':\n{}\n{}\n".format(ret, out, err))

def _run_experiment(args):
    """Converts the RDF data into a machine-learnable dataset.

    :param args.target: target experiment.
    :param args.dst: destination directory for the results.
    :param args.endpoint: Virtuoso SPARQL endpoint.
    :param args.default_graph: URI of the default graph.
    """
    ALL_TARGETS = {
        "yip" : experiments.YipExperiment,
        "sgd" : experiments.SGDExperiment,
    }

    go_aspects = None
    if args.go_aspects is not None:
        go_aspects = args.go_aspects.split(",")

    endpoint = ocelot.Endpoint(args.endpoint, args.default_graph)

    for id_ in ocelot.validate(ALL_TARGETS.keys(), args.targets):
        Experiment = ALL_TARGETS[id_]
        experiment = Experiment(args.src, args.dst, endpoint=endpoint,
                                go_aspects=go_aspects,
                                max_go_depth=args.max_go_depth,
                                min_go_annot=args.min_go_annot,
                                rng=args.seed)
        experiment.run(force=args.force_update)

def main():
    import argparse as ap

    COMMANDS = {
        "make-rdf": _make_rdf,
        "start-virtuoso": _start_virtuoso,
        "upload-rdf": _upload_rdf,
        "clear-graph": _clear_graph,
        "run-experiment": _run_experiment,
    }

    parser = ap.ArgumentParser(description="Ocelot, the playful feline (also SBR for proteomics).")
    parser.add_argument("command", type=str, default=None,
                        help="What to do, any of {}".format(COMMANDS.keys()))
    parser.add_argument("-t", "--targets", nargs="+",
                        help="command-specific targets (default: all)")
    parser.add_argument("-s", "--src", type=str,
                        help="path to the source data")
    parser.add_argument("-d", "--dst", type=str,
                        help="path to the results destination")
    parser.add_argument("-e", "--endpoint", type=str,
                        default = ocelot.config["rdf"][u"endpoint"],
                        help="URI of SPARQL endpoint to use (default: value in ocelot.ini)")
    parser.add_argument("-g", "--default-graph", type=str,
                        default = ocelot.config["rdf"][u"default_graph"],
                        help="URI of the default graph (default: value in ocelot.ini)")
    parser.add_argument("-f", "--force-update", action="store_true",
                        help="ignores cached data")
    parser.add_argument("--debug", action="store_true",
                        help="[all] enable debugging output")
    parser.add_argument("--seed", type=int, default=None,
                        help="[all] seed for the RNG.")
    parser.add_argument("--go-aspects", type=str, default=None,
                        help="[run-experiment] restrict GO annotations/predictions to these aspects (default: all of them)")
    parser.add_argument("--max-go-depth", type=int, default=None,
                        help="[run-experiment] restrict GO annotations/predictions to terms with at most this depth (default: None)")
    parser.add_argument("--min-go-annot", type=int, default=None,
                        help="[run-experiment] restrict GO annotations/predictions to term with at least this many annotations (default: None)")

    args = parser.parse_args()

    command = COMMANDS.get(args.command)
    if not command:
        raise ValueError("invalid command '{}'".format)
    command(args)

if __name__ == "__main__":
    main()
