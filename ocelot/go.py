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
    alt_ids : list
        List of alternative IDs.
    name : str
        Short description of the term.
    namespace : str
        Namespace of the term, one of ``["biological_process",
        "cellular_component", "molecular_function"]``.
    is_obsolete : bool
        Whether the term is marked as obsolete.
    level : set
        Distance from the root, ``-1`` if unknown.
    parents : list
        List of (parent term ID, relation) tuples.
    children : list
        List of (child term ID, relation) tuples.
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
        self.parents = []
        self.children = []
        self.proteins = set()

    def __str__(self):
        return "{}\tlevel {}\t{} [{}] {}" \
                    .format(self.id, self.level, self.name, self.namespace,
                            "**OBSOLETE**" if self.is_obsolete else "")

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def get_parents(self, dag, relations=set(["is_a"])):
        """Returns the parents of the node."""
        parents = set()
        for parent_id, relation in self.parents:
            if not relation in relations:
                continue
            parents.add((dag._id_to_term[parent_id], relation))
        return parents

    def get_children(self, dag, relations=set(["is_a"])):
        """Returns the children of the term."""
        children = set()
        for child_id, relation in self.children:
            if not relation in relations:
                continue
            children.add((dag._id_to_term[child_id], relation))
        return children

    def get_ancestors(self, dag, relations=set(["is_a"])):
        """Returns the ancestors of the node.

        :param dag: DAG.
        :param relations:
        :returns: ``set`` of term IDs.
        """
        ancestors = set()
        for parent_id, relation in self.parents:
            if not relation in relations:
                continue
            ancestors.add(dag._id_to_term[parent_id].id)
            ancestors |= dag._id_to_term[parent_id].get_ancestors(dag, relations)
        return ancestors

    def get_descendants(self, dag, relations=set(["is_a"])):
        """Returns the descendants of a node."""
        descendants = set()
        for child_id, relation in self.parents:
            if not relation in relations:
                continue
            descendants.add(dag._id_to_term[child_id].id)
            descendants |= dag._id_to_term[child_id].get_descendants(dag, relations)
        return descendants

class BinTerm(GOTerm):
    """A placeholder GO term."""

    def __init__(self, parent_id, parent_term, proteins=set()):
        super(BinTerm, self).__init__()

        self.id = "BN:{}".format(parent_id[parent_id.index(":")+1:])
        self.name = "bin term for {}".format(parent_id)
        self.namespace = parent_term.namespace
        self.is_obsolete = False
        self.level = parent_term.level + 1
        self.parents = [(parent_id, "bin_for")]
        self.children = []
        self.proteins = proteins

class _GOReader(object):
    """Parse a GO OBO file.

    The latest version of the OBO files can be downloaded from::

        http://purl.obolibrary.org/obo/go/

    Parameters
    ----------
    path : str
        Path to the OBO file.
    keep_obsolete : bool, defaults to ``True``
        Whether to keep obsolete nodes.
    """
    _TYPEDEF_TAG    = "[Typedef]"
    _TERM_TAG       = "[Term]"

    def __init__(self, path, keep_obsolete = True):
        self._handle = open(path)
        self._keep_obsolete = keep_obsolete

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
            while not done:
                term = self.next()
                done = not term.is_obsolete
            yield term

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
                term.parents.append((value.split()[0], "is_a"))
            elif key == "relationship":
                term.parents.append((value.split()[1], value.split()[0]))
            elif key == "is_obsolete":
                term.is_obsolete = {"true":True, "false":False}[value]
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
    def __init__(self, path=None, keep_obsolete=False, include_alt_ids=False):
        self._id_to_term = {}
        self._path = path
        if not path:
            return
        for term in _GOReader(path, keep_obsolete = keep_obsolete):
            assert not term.id in self._id_to_term
            self._id_to_term[term.id] = term
            if include_alt_ids:
                for alt_id in term.alt_ids:
                    assert not alt_id in self._id_to_term
                    self._id_to_term[alt_id] = term
        self._populate_terms()

    def _set_term_depth(self, term):
        if term.level < 0:
            if not term.parents:
                term.level = 0
            else:
                parent_depths = [self._set_term_depth(self._id_to_term[parent])
                                 for parent, relation in term.parents
                                 if parent in self._id_to_term and relation in ("is_a", "bin_for")]
                if not len(parent_depths):
                    term.level = -1
                else:
                    term.level = min(parent_depths) + 1
        return term.level

    def _populate_terms(self):
        for term in self._id_to_term.itervalues():
            for parent, relation in term.parents:
                self._id_to_term[parent].children.append((term.id, relation))
            if term.level < 0:
                self._set_term_depth(term)

    def __repr__(self):
        return "GODag('{}')".format(self._path)

    def __str__(self):
        return "\n".join("{}: {}".format(id_, term)
                         for id_, term in self._id_to_term.items())

    def __getitem__(self, term_id):
        """Returns the ``GOTerm`` associated to a term ``id``."""
        if not term_id in self._id_to_term:
            return
        return self._id_to_term[term_id]

    def paths_to_root(self, term_id):
        """Returns all possible paths to the root node.

        Each path includes the term given. The order of the path is
        top -> bottom, i.e. it starts with the root and ends with the
        given term (inclusively).

        :param term_id: id of the GO term (e.g., ``'GO:0003682'``).
        :returns: a list of GO terms.
        """
        if term_id not in self._id_to_term:
            return

        def _paths_to_root(term):
            if term.level == 0:
                return [[term]]
            paths = []
            for parent_id, relation in term.parents:
                paths_above = _paths_to_root(self._id_to_term[parent_id])
                for path_above in paths_above:
                    paths.append([ term ] + path_above)
            return paths

        return _paths_to_root(self._id_to_term[term_id])

    def get_valid_term_ids(self, include_bins=False):
        for id_ in self._id_to_term.iterkeys():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_

    def get_valid_terms(self, include_bins=False):
        for id_, term in self._id_to_term.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield term

    def get_valid_ids_terms(self, include_bins=False):
        for id_, term in self._id_to_term.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_, term

    def get_interaspect_links(self):
        """Generates the inter-aspect relations.

        :returns: triples of the form (term, relation, parent).
        """
        for term in self._id_to_term.itervalues():
            for parent_id, relation in term.parents:
                parent = self._id_to_term[parent_id]
                if term.namespace != parent.namespace:
                    yield term, relation, parent

    def _check_annotations(self, p_to_term_ids):
        num_annotations_1 = sum(len(term_ids)
                                for term_ids in p_to_term_ids.itervalues())
        num_annotations_2 = sum(len(term.proteins)
                                for term in self._id_to_term.itervalues())
        return num_annotations_1 == num_annotations_2, "{} != {}".format(num_annotations_1, num_annotations_2)

    def _check_annotation_hierarchy(self):
        for term in self._id_to_term.itervalues():

            # Sanity check: check that a term has at most as many annotations
            # as the number of annotations in all its parents (note that a term
            # may have more than one parent, hence the sum).
            num_parent_annot = 0
            for parent_id, relation in term.parents:
                parent = self._id_to_term[parent_id]
                if relation == "is_a":
                    num_parent_annot += len(parent.proteins)
            assert len(term.parents) == 0 or num_parent_annot >= len(term.proteins), \
                "failed sanity check: {} -> parents {}".format(term, term.parents)

            # Sanity check: check that a term has at least as many annotations
            # as each of its children
            for child_id, relation in term.children:
                child = self._id_to_term[child_id]
                if relation == "is_a":
                    assert len(term.proteins) >= len(child.proteins), \
                        "failed sanity check: {} -> children {}".format(term, term.children)

    def _check(self, p_to_term_ids):
        assert sum(len(term.proteins) for term in self._id_to_term.itervalues())
        assert self._check_annotations(p_to_term_ids)
        self._check_annotation_hierarchy()

    def _propagate(self, p_to_term_ids):
        """Propagate annotations to the root."""
        propagated_p_to_term_ids = {}
        for p, term_ids in p_to_term_ids.iteritems():
            propagated_term_ids = set()
            for term_id in term_ids:
                if term_id == "":
                    print "Warning: protein '{}' is annotated with an empty term, skipping" \
                            .format(p)
                    continue
                paths = self.paths_to_root(term_id)
                if paths is None:
                    print "Warning: protein '{}' is annotated with an unknown term ID '{}', skipping" \
                            .format(p, term_id)
                    continue
                for path in paths:
                    # Do the actual propagation
                    for term in path:
                        term.proteins.add(p)
                        print "propagating", term, "with", p
                    # Take note of the propagated terms
                    propagated_term_ids.update(term.id for term in path)
            # Update the protein ID->term IDs map
            propagated_p_to_term_ids[p] = propagated_term_ids

        self._check(propagated_p_to_term_ids)

        return propagated_p_to_term_ids

    def annotate(self, p_to_term_ids, propagate = False):
        """Annotates the GO DAG with protein annotations.

        :param p_to_term_ids: map from protein IDs to GO Term IDs.
        :param propagate: whether to propagate the annotations to the root.
        :returns: the annotated GO DAG.
        """
        for p, term_ids in p_to_term_ids.iteritems():
            for term_id in term_ids:
                assert term_id in self._id_to_term, "unknown term ID '{}'".format(term_id)
                self._id_to_term[term_id].proteins.add(p)
                print "annotating", term_id, "with", p

        self._check(p_to_term_ids)

        if propagate:
            p_to_term_ids = self._propagate(p_to_term_ids)

        return p_to_term_ids

    def get_proteins_by_aspect(self):
        """Returns the proteins that only have some aspects annotated."""
        aspect_to_ps = defaultdict(set)
        for term_id in self._id_to_term:
            term = self._id_to_term[term_id]
            aspect_to_ps[term.namespace].update(term.proteins)
        return aspect_to_ps

    def generate_terms_by_level(self):
        """Generates all terms according to their recorder level."""
        level_to_terms = defaultdict(set)
        for term in self._id_to_term.itervalues():
            level_to_terms[term.level].add(term)
        levels = sorted(level_to_terms.keys())
        for level in levels:
            for term in level_to_terms[level]:
                yield term

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

        self._check_annotation_hierarchy()

        # Do a per-level traversal of the dag, and mark the terms that satisfy
        # all constraints.
        processed, terms_to_keep = set(), set()
        for term in self.generate_terms_by_level():
            if term in processed:
                continue
            processed.add(term)

            if not max_depth is None and term.level > max_depth:
                print "discarding '{}', too deep ({} > {})".format(term, term.level, max_depth)
                continue
            if not term.namespace in aspects:
                print "discarding '{}', not in namespace".format(term)
                continue
            if not min_annot is None and len(term.proteins) < min_annot:
                print "discarding '{}', too few annotations ({} < {})".format(term, len(term.proteins), min_annot)
                continue

            terms_to_keep.add(term)

        print "found {} terms to keep".format(len(terms_to_keep))

        # Ancestores of terms to be kept must also be kept
        expanded_terms_to_keep = set()
        for term in terms_to_keep:
            expanded_terms_to_keep.update(self._id_to_term[id_] for id_ in term.get_ancestors(self))
        terms_to_keep = terms_to_keep | expanded_terms_to_keep

        print "expanded to {} terms to keep".format(len(terms_to_keep))

        # Sanity check
        for term in terms_to_keep:
            if not max_depth is None and term.level > max_depth:
                print "Warning: keeping '{}' even if too deep ({} > {})".format(term, term.level, max_depth)
            if not term.namespace in aspects:
                print "Warning: keeping '{}' even if not in allowed namespaces".format(term)
            if not min_annot is None and len(term.proteins) < min_annot:
                print "Warning: keeping '{}' even if too few annotations ({} < {})".format(term, len(term.proteins), min_annot)

        # Create a new id->term dictionary for terms to be kept only
        filtered_id_to_term = {term.id: deepcopy(term) for term in terms_to_keep}

        # Remove all relations to terms to be discarded
        for term in filtered_id_to_term.itervalues():
            term.parents = [(id_, rel) for id_, rel in term.parents
                            if id_ in filtered_id_to_term]
            term.children = [(id_, rel) for id_, rel in term.children
                             if id_ in filtered_id_to_term]
            term.level = -1

        print "len(filtered_id_to_term) = {}".format(len(filtered_id_to_term))

        # Sanity check
        for term in filtered_id_to_term.itervalues():
            assert all(id_ in filtered_id_to_term for id_, _ in term.parents), \
                   "invalid parents for {}:\n{}".format(term, term.parents)
            assert all(id_ in filtered_id_to_term for id_, _ in term.children), \
                   "invalid children for {}:\n{}".format(term, term.children)

        # Add bin nodes as children of any term which has some retained and
        # some removed children; the proteins of the removed chilren are
        # assigned to the new bin term
        bin_terms = set()
        for id_ in filtered_id_to_term:
            removed_children = set(self._id_to_term[id_].children) - set(filtered_id_to_term[id_].children)

            removed_proteins = set()
            for child_id, _ in removed_children:
                removed_proteins.update(self._id_to_term[child_id].proteins)

            if len(removed_proteins) > (min_annot if min_annot is None else 0):
                print "adding bin term to {}".format(id_)
                bin_term = BinTerm(id_, filtered_id_to_term[id_],
                                   proteins=removed_proteins)
                assert not bin_term.id in bin_terms, \
                       "multiple insertions of bin term {}".format(bin_term.id)
                bin_terms.add(bin_term)

        for bin_term in bin_terms:
            filtered_id_to_term[bin_term.id] = bin_term

        print "added {} bin terms".format(len(bin_terms))

        # Sanity check
        for term in filtered_id_to_term.itervalues():
            assert all(id_ in filtered_id_to_term for id_, _ in term.parents), \
                   "invalid parents for {}:\n{}".format(term, term.parents)
            assert all(id_ in filtered_id_to_term for id_, _ in term.children), \
                   "invalid children for {}:\n{}".format(term, term.children)

        # Make sure that no unwanted term ended up in the filtered list
        ids_to_keep = set(term.id for term in terms_to_keep)
        kept_ids_minus_bins = set(id_ for id_ in filtered_id_to_term if not id_.startswith("BN:"))
        assert ids_to_keep == kept_ids_minus_bins

        # Compute the new protein->term map
        filtered_p_to_term_ids = defaultdict(set)
        for id_, term in filtered_id_to_term.iteritems():
            for p in term.proteins:
                filtered_p_to_term_ids[p].add(id_)

        # Make sure that all original proteins are still mapped to at least
        # one kept term
        assert all(p in filtered_p_to_term_ids and len(filtered_p_to_term_ids[p]) for p in ps)

        # Replace the current id->term map
        self._id_to_term = filtered_id_to_term

        # Recompute the term depth
        for term in self._id_to_term.itervalues():
            self._set_term_depth(term)

        self._check_annotation_hierarchy()

        return filtered_p_to_term_ids

    def draw(self, path, p_to_term_ids, fontname="Verdana"):
        import pydot

        NAMESPACE_TO_NS = {
            "biological_process": "bp",
            "cellular_component": "cc",
            "molecular_function": "mf",
        }

        graph = pydot.Dot(graph_type="digraph", fontname=fontname)

        term_to_node = {}
        for _, term_ids in p_to_term_ids.items():
            for term_id in term_ids:
                term = self._id_to_term[term_id]
                name = term.id.replace(":", "") + "_" + NAMESPACE_TO_NS[term.namespace]
                node = pydot.Node(name, fontname=fontname)
                term_to_node[term] = node
                graph.add_node(node)

        for _, term_ids in p_to_term_ids.items():
            for term_id in term_ids:
                term = self._id_to_term[term_id]
                for parent_id, relation in term.parents:
                    parent = self._id_to_term[parent_id]
                    graph.add_edge(pydot.Edge(term_to_node[term], term_to_node[parent], fontname=fontname))

        graph.write_png(path)

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
