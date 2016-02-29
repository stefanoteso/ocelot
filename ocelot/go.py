# -*- coding: utf-8 -*-

import sys, re, copy
import numpy as np
from collections import defaultdict
from ocelot.services import _cls
from copy import deepcopy

class GOTerm(object):
    """A GO Term.

    Attributes
    ----------
    id : str
        ID of the term, e.g. ``"GO:0048308'"``
    alt_ids : list of str
        List of alternative IDs.
    name : str
        Short description of the term.
    namespace : str
        Namespace of the term, one of ``["biological_process",
        "cellular_component", "molecular_function"]``.
    is_obsolete : bool
        Whether the term is marked as obsolete.
    level : int
        Distance from the root, ``-1`` if unknown.
    proteins : set
        Collection of proteins annotated with the GO term.
    """
    def __init__(self):
        self.id = ""
        self.alt_ids = []
        self.name = ""
        self.namespace = ""
        self.is_obsolete = False
        self.level = -1
        self.proteins = set()
        self._parents = []
        self._children = []

    def __str__(self):
        return "{}\tlevel {}\t{} [{}] {}" \
                    .format(self.id, self.level, self.name, self.namespace,
                            "**OBSOLETE**" if self.is_obsolete else "")

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def get_parents(self, dag, relations=frozenset(["is_a"])):
        """Computes all parents of the current node wrt the DAG.

        It only returns parents according to the specified relations, or all
        of them if relations is None.

        :returns: a set of (term, relation string) pairs.
        """
        return set((dag._id_to_term[id_], relation)
                   for id_, relation in self._parents
                   if relations is None or relation in relations)

    def get_children(self, dag, relations=frozenset(["is_a"])):
        """Computes all children of the current node wrt the DAG.

        It only returns children according to the specified relations, or all
        of them if relations is None.

        :returns: a set of (term, relation string) pairs.
        """
        return set((dag._id_to_term[id_], relation)
                   for id_, relation in self._children
                   if relations is None or relation in relations)

    def get_paths_to_root(self, dag, relations=frozenset(["is_a"])):
        """Computes all paths from the term to the root.

        :returns: a list of paths.
        """
        def recurse(term, dag, relations):
            parents = term.get_parents(dag, relations=relations)
            if len(parents) == 0:
                return [[term]]
            paths = []
            for parent, _ in parents:
                paths_above = recurse(parent, dag, relations)
                for path in paths_above:
                    paths.append(path + [term])
            return paths

        return recurse(self, dag, relations)

    def update_level(self, dag, relations=frozenset(["is_a"])):
        """Computes the term level and updates it.

        It also updates the level of all parents.

        :returns: the level.
        """
        parents = self.get_parents(dag, relations=relations)
        if not len(parents):
            self.level = 0
        else:
            parent_levels = [parent.update_level(dag)
                             for parent, _ in self.get_parents(dag, relations=relations)]
            self.level = min(parent_levels) + 1
        return self.level

class BinTerm(GOTerm):
    """A placeholder GO term."""

    def __init__(self, parent_id, parent, proteins):
        super(BinTerm, self).__init__()

        self.id = "BN:{}".format(parent_id.split(":", 1)[-1])
        self.name = "bin term for {}".format(parent_id)
        self.namespace = parent.namespace
        self.is_obsolete = parent.is_obsolete
        self.level = parent.level + 1
        self._parents = [(parent_id, "bin_for")]
        self._children = tuple()
        self.proteins = proteins

class _GOReader(object):
    """Parse a GO OBO file.

    The latest version of the OBO files can be downloaded from::

        http://purl.obolibrary.org/obo/go/

    Parameters
    ----------
    path : str
        Path to the OBO file.
    """
    _TYPEDEF_TAG    = "[Typedef]"
    _TERM_TAG       = "[Term]"

    def __init__(self, path):
        self._handle = open(path)

    @staticmethod
    def _read_until(handle, prefix):
        """Read each line until a prefix is found; then puts the prefix back."""
        from exceptions import EOFError
        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(prefix):
                handle.seek(pos)
                return
        raise EOFError("'{}' prefix not found".format(prefix))

    def __iter__(self):
        line = self._handle.readline()
        if not line.startswith(self._TERM_TAG):
            self._read_until(self._handle, self._TERM_TAG)
        while True:
            done = False
            yield self.next()

    def next(self):
        line = self._handle.readline()
        if not line or line.startswith(self._TYPEDEF_TAG):
            raise StopIteration

        # Read until the next tag and save everything in between
        lines = []
        while True:
            pos = self._handle.tell()
            line = self._handle.readline()
            if not line or \
                (line.startswith(self._TYPEDEF_TAG) or \
                 line.startswith(self._TERM_TAG)):
                self._handle.seek(pos)
                break
            lines.append(line)

        # XXX the full set of keys is: alt_id, auto-generated-by, comment,
        # consider, created_by, creation_date, data-version, date, def,
        # default-namespace, disjoint_from, expand_assertion_to,
        # format-version, holds_over_chain, id, intersection_of, inverse_of,
        # is_a, is_class_level, is_metadata_tag, is_obsolete, is_transitive,
        # name, namespace, ontology, property_value, relationship, remark,
        # replaced_by, saved-by, subset, subsetdef, synonym, synonymtypedef,
        # transitive_over, xref. We do not handle nearly enough of them.
        #
        # XXX the full set of relations is: happens_during, has_part,
        # negatively_regulates, occurs_in, part_of, positively_regulates,
        # regulates.

        term = GOTerm()
        for line in lines:
            match = re.match("^([^ ]+): (.+)$", line)
            if not match:
                continue
            assert len(match.groups()) == 2
            key, value = match.groups()
            if key == "id":
                term.id = value
            elif key == "alt_id":
                term.alt_ids.append(value)
            elif key == "name":
                term.name = value
            elif key == "namespace":
                term.namespace = value
            elif key == "is_a":
                term._parents.append((value.split()[0], "is_a"))
            elif key == "relationship":
                term._parents.append((value.split()[1], value.split()[0]))
            elif key == "is_obsolete":
                term.is_obsolete = {"true": True, "false": False}[value]
        return term

class GODag(object):
    """The GO DAG.

    Essentially a map from term IDs and GOTerm objects.

    Parameters
    ----------
    path : str
        Path to the OBO file.
    keep obsolete : bool, defaults to ``False``
        Whether to keep terms marked as obsolete.
    include_alt_ids : bool, defaults to ``False``
        Whether to include alternative term IDs in the dictionary.
    """
    def __init__(self, path, keep_obsolete=False, include_alt_ids=False):
        self._id_to_term = {}
        self._path = path

        # Read the DAG from the OBO file
        for term in _GOReader(path):
            assert not term.id in self._id_to_term

            if term.is_obsolete and not keep_obsolete:
                print "Warning: skipping obsolete GO term '{}'".format(term.id)
                continue

            self._id_to_term[term.id] = term
            if include_alt_ids:
                for alt_id in term.alt_ids:
                    assert not alt_id in self._id_to_term
                    self._id_to_term[alt_id] = term

        # Compute the children nodes
        for term in self._id_to_term.itervalues():
            for parent, relation in term.get_parents(self):
                parent._children.append((term.id, relation))

        # Update the term levels
        for term in self._id_to_term.itervalues():
            term.update_level(self)

    def __repr__(self):
        return "GODag('{}')".format(self._path)

    def get_id_to_term(self):
        return self._id_to_term

    def get_p_to_term_ids(self, id_to_term=None):
        """Computes the protein -> term ID map."""
        if id_to_term is None:
            id_to_term = self._id_to_term

        p_to_term_ids = defaultdict(set)
        for term in id_to_term.itervalues():
            for p in term.proteins:
                p_to_term_ids[p].add(term.id)
        return dict(p_to_term_ids)

    def get_terms_by_level(self):
        level_to_terms = defaultdict(set)
        for term in self._id_to_term.itervalues():
            level_to_terms[term.level].add(term)
        for level in sorted(level_to_terms.keys()):
            for term in level_to_terms[level]:
                yield term

    def draw(self, path, fontname="Verdana"):
        """Draws the annotated DAG to a PNG file."""
        import pydot

        NAMESPACE_TO_NS = {
            "biological_process": "bp",
            "cellular_component": "cc",
            "molecular_function": "mf",
        }

        p_to_term_ids = self.get_p_to_term_ids()

        graph = pydot.Dot(graph_type="digraph", fontname=fontname)

        term_to_node = {}
        for _, term_ids in p_to_term_ids.items():
            for term_id in term_ids:
                term = self._id_to_term[term_id]
                name = NAMESPACE_TO_NS[term.namespace] + "_" + term.id.replace(":", "") + "_" + str(len(term.proteins))
                node = pydot.Node(name, fontname=fontname)
                term_to_node[term] = node
                graph.add_node(node)

        for _, term_ids in p_to_term_ids.items():
            for term_id in term_ids:
                term = self._id_to_term[term_id]
                for parent, _ in term.get_parents(self):
                    graph.add_edge(pydot.Edge(term_to_node[term], term_to_node[parent], fontname=fontname))

        graph.write_png(path)

    def _check_term_parent_annotations(self):
        for term in self._id_to_term.itervalues():

            # All terms must have at *most* as many annotations as any of their
            # parents
            parents = [p for p, r in term.get_parents(self)]
            for parent in parents:
                assert len(term.proteins) <= len(parent.proteins), \
                    "term '{}' has more proteins ({}) than its parent '{}' ({})" \
                        .format(term, len(term.proteins), parent, len(parent.proteins))

            # All terms must have at *least* as many annotations as any of
            # their children
            children = [c for c, r in term.get_children(self)]
            for child in children:
                assert len(child.proteins) <= len(term.proteins), \
                    "term '{}' has less proteins ({}) than its child '{}' ({})" \
                        .format(term, len(term.proteins), child, len(child.proteins))

    def _check_dag(self, id_to_term, aspects, max_depth, min_annot):
        num_roots = 0
        for id_, term in id_to_term.iteritems():

            # All terms must respect the requirements
            assert term.namespace in aspects
            assert term.level >= 0
            assert max_depth is None or term.level <= max_depth
            assert min_annot is None or len(term.proteins) >= min_annot

            # All parents (over any relation) of a term must also be in the
            # id_to_term map
            assert all(parent.id in id_to_term or parent.namespace != term.namespace
                       for parent, _ in term.get_parents(self, relations=None))

            # All children (over any relation) of a term must also be in the
            # id_to_term map
            assert all(child.id in id_to_term or child.namespace != term.namespace
                       for child, _ in term.get_children(self, relations=None))

            # XXX this also matches obsolete terms!
            if term.level == 0:
                num_roots += 1

        # There must be exactly one root per aspect
        assert num_roots == len(aspects), "roots {} mismatch with aspects {}"\
            .format(num_roots, len(aspects))

    def annotate(self, p_to_term_ids, propagate=False):
        """Annotates the GO DAG with protein annotations.

        Warning: unknown annotations (e.g. referring to obsolete GO terms, if
        ``keep_obsolete`` is ``False``) are not retained. Consequent calls to
        get_p_to_term_ids() will not includes those annotations.

        Warning: annotations are only propagated through ``is_a`` relations; it
        is not meaningful to propagate over the other relations. IOW, for the
        purpose of this method, we treat the GO DAG as nothing more than a
        set/subset hiearchy.

        Parameters
        ----------
        p_to_term_ids : dict
            Map from protein ID to GO term IDs.
        propagate : bool, defaults to ``False``
            Whether to propagate the annotations to the root.
        """
        for p, term_ids in p_to_term_ids.iteritems():
            for term_id in term_ids:
                if not term_id in self._id_to_term:
                    print "Warning: protein '{}' is annotated with an unknown term ID '{}'" \
                            .format(p, term_id)
                else:
                    self._id_to_term[term_id].proteins.add(p)

        if not propagate:
            return

        for term in self._id_to_term.values():
            if not len(term.proteins):
                continue

            paths = term.get_paths_to_root(self)
            assert len(paths), "term '{}' has no paths".format(term.id)

            for path in paths:
                for ancestor in path:
                    ancestor.proteins.update(term.proteins)

        self._check_term_parent_annotations()

    def preprocess(self, ps, aspects=None, max_depth=None, min_annot=None):
        """Processes the DAG by removing unwanted terms and adding bin terms.

        This method does the following:

        * Removes all terms in unwanted namespaces.
        * Removes all terms without enough annotations.
        * Removes all terms that are too deep.
        * Adds *bin* terms.

        The protein-to-function map is adjusted accordingly.

        Parameters
        ----------
        ps : list
            List of protein IDs a strings
        aspects : list, optional, defaults to ``None``
            Terms not in these aspects are discarded.
        max_depth : int, optional, defaults to ``None``
            Terms deeper than this are discarded.
        min_annot : int, optional, defaults to ``None``
            Terms with fewer annotations than this are discarded.
        """
        ALL_ASPECTS = [
            "biological_process",
            "cellular_component",
            "molecular_function"
        ]

        if aspects is None:
            aspects = ALL_ASPECTS
        elif any(aspect not in ALL_ASPECTS for aspect in aspects):
            raise ValueError("invalid GO aspect.")

        # Store all terms that we wish to keep
        ids_to_keep = set()
        for term in self._id_to_term.itervalues():
            if term.level < 0:
                print "discarding '{}', invalid level".format(term.id)
                continue
            if not term.namespace in aspects:
                print "discarding '{}', not in namespace".format(term.id)
                continue
            if not max_depth is None and term.level > max_depth:
                print "discarding '{}', too deep ({} > {})".format(term.id, term.level, max_depth)
                continue
            if not min_annot is None and len(term.proteins) < min_annot:
                print "discarding '{}', too few annotations ({} < {})".format(term.id, len(term.proteins), min_annot)
                continue
            ids_to_keep.add(term.id)

        # Create a new id->term dictionary for terms to be kept only and prune
        # their topology. Even though the above pruning is only based on the
        # ``is_a`` relation, here we preserve also the other relations.
        filtered_id_to_term = {}
        for id_ in ids_to_keep:
            term = deepcopy(self._id_to_term[id_])
            term._parents = [(parent_id, relation)
                             for parent_id, relation in term._parents
                             if parent_id in ids_to_keep]
            term._children = [(child_id, relation)
                              for child_id, relation in term._children
                              if child_id in ids_to_keep]
            filtered_id_to_term[id_] = term

        self._check_dag(filtered_id_to_term, aspects, max_depth, min_annot)

        # Add the bin nodes. A bin is added whenever at least one child of a
        # retained term has been removed and and at least one child has been
        # kept. No bin node is added if the removed child had an insufficent
        # number of proteins, nor when the maximum term depth has been reached.
        bin_terms, bin_ids = set(), set()
        for id_ in filtered_id_to_term.iterkeys():
            old_term = self._id_to_term[id_]
            old_children_ids = set(c.id for c, r in old_term.get_children(self))

            new_term = filtered_id_to_term[id_]
            new_children_ids = set(c.id for c, r in new_term.get_children(self))

            assert len(new_children_ids - old_children_ids) == 0

            # Collect annotations in removed terms
            removed_proteins = set()
            for child_id in old_children_ids - new_children_ids:
                removed_proteins.update(self._id_to_term[child_id].proteins)

            # Check if we must add a bin term
            if len(removed_proteins) < (1 if min_annot is None else min_annot):
                continue
            if max_depth is not None and new_term.level == max_depth:
                continue

            # Add the bin term
            bin_term = BinTerm(id_, new_term, removed_proteins)
            assert not bin_term.id in bin_ids, \
                "multiple insertion of the same bin term '{}'".format(bin_term.id)

            print "Adding bin term '{}' (with {} proteins) to '{}'" \
                    .format(bin_term.id, len(bin_term.proteins), id_)
            bin_terms.add(bin_term)
            bin_ids.add(bin_term.id)

        # Update the id to term map
        for bin_term in bin_terms:
            filtered_id_to_term[bin_term.id] = bin_term

        print "Added {} bin terms".format(len(bin_terms))

        self._check_dag(filtered_id_to_term, aspects, max_depth, min_annot)

        # Replace the term ID -> term map
        self._id_to_term = filtered_id_to_term


class _TestDAG(object):

    def _load_and_fill(self, p_to_term_ids, propagate):
        dag = GODag("data/GO/go-basic.obo")
        return dag, dag.annotate(p_to_term_ids, propagate=propagate)

    def _all_term_ids(self, p_to_term_ids):
        all_term_ids = set()
        for term_ids in p_to_term_ids.values():
            all_term_ids.update(term_ids)
        return all_term_ids

    def test_propagation_roots(self):
        ANNOTATIONS = {
            "BP": ["GO:0008150"], # biological process
            "MF": ["GO:0003674"], # molecular function
            "CC": ["GO:0005575"], # cellular component
        }

        for propagate in (True, False):
            dag, p_to_term_ids = self._load_and_fill(ANNOTATIONS, propagate)

            assert len(p_to_term_ids) == 3
            assert sum(len(term.proteins) for _, term in dag._id_to_term.items()) == 3

    def test_propagation_first_level(self):
        ANNOTATIONS = {
            "p0": ["GO:0009987"], # cellular process
            "p1": ["GO:0032502"], # developmental process
            "p2": ["GO:0002376"], # immune system process
            "p3": ["GO:0008152"], # metabolic process
            "p4": ["GO:0051704"], # multi-organism process
            "p5": ["GO:0032501"], # multicellular organismal process
            "p6": ["GO:0022414"], # reproductive process
            "p7": ["GO:0048511"], # rhythmic process
            "p8": ["GO:0044699"], # single-organism process
        }

        _, p_to_term_ids = self._load_and_fill(ANNOTATIONS, False)
        assert len(self._all_term_ids(p_to_term_ids)) == 9

        _, p_to_term_ids = self._load_and_fill(ANNOTATIONS, True)
        assert len(self._all_term_ids(p_to_term_ids)) == 10
