# -*- coding: utf-8 -*-

import sys, re, copy
from collections import defaultdict
from ocelot.services import _cls

class _GOReader(object):
    """Parse a GO OBO file.

    The latest version of the OBO files can be downloaded from::

        http://purl.obolibrary.org/obo/go/

    :param path: path to the OBO file.
    :param keep_obsolete: whether to keep obsolete nodes (default: ``True``).

    Example usage::

        >>> reader = GOReader()
        >>> for term in reader:
                print term
    .. todo::

        Handle the remaining term information.
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

class GOTerm(object):
    """A GO Term."""
    def __init__(self):
        self.id = ""
        """ID of the term, e.g. ``"GO:0048308'"``"""
        self.alt_ids = []
        """List of alternative IDs."""
        self.name = ""
        """Short description of the term."""
        self.namespace = ""
        """Namespace of the term, one of ``["biological_process", "cellular_component", "molecular_function"]``."""
        self.parents = []
        """List of (parent term ID, relation) tuples."""
        self.children = []
        """List of (child term ID, relation) tuples."""
        self.is_obsolete = False
        """Whether the term is marked as obsolete."""
        self.level = -1
        """Distance from the root, ``-1`` if unknown."""
        self.proteins = set()

    def __str__(self):
        return "{}\tlevel {}\t{} [{}] {}" \
                    .format(self.id, self.level, self.name, self.namespace,
                            "**OBSOLETE**" if self.is_obsolete else "")

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def get_name(self):
        bound = 15
        ret = self.name[0].upper() + self.name[1:].replace('_', ' ')
        if len(ret) > bound:
            ret = ret[:bound] + ret[bound:].replace(' ', '\n', 1)
        return ret

    def has_ancestor(self, dag, term_id, relations=set(["is_a"])):
        """Checks whether ``term_id`` is an ancestor of this term."""
        for parent_id, relation in self.parents:
            if not relation in relations:
                continue
            if dag._terms[parent_id].id == term_id or \
               dag._terms[parent_id].has_parent(term_id):
                return True
        return False

    def has_descendant(self, dag, term, relations=set(["is_a"])):
        """Checks whether ``term_id`` is a descendant of this term."""
        for child_id, relation in self.children:
            if not relation in relations:
                continue
            if dag._terms[child_id].id == term_id or \
               dag._terms[child_id].has_child(term_id):
                return True
        return False

    def get_parents(self, dag, relations=set(["is_a"])):
        """Returns the parents of a node."""
        parents = set()
        for parent_id, relation in self.parents:
            if not relation in relations:
                continue
            parents.add((dag[parent_id], relation))
        return parents

    def get_children(self, dag, relations=set(["is_a"])):
        """Returns the children of a node."""
        children = set()
        for child_id, relation in self.children:
            if not relation in relations:
                continue
            children.add((dag[child_id], relation))
        return children

    def get_ancestors(self, dag, relations=set(["is_a"])):
        """Returns the ancestors of a node."""
        ancestors = set()
        for parent_id, relation in self.parents:
            if not relation in relations:
                continue
            ancestors.add(dag._terms[parent_id].id)
            ancestors |= dag._terms[parent_id].get_ancestors(dag, relations)
        return ancestors

    def get_descendants(self, dag, relations=set(["is_a"])):
        """Returns the descendants of a node."""
        descendants = set()
        for child_id, relation in self.parents:
            if not relation in relations:
                continue
            descendants.add(dag._terms[child_id].id)
            descendants |= dag._terms[child_id].get_descendants(dag, relations)
        return descendants

class GODag(object):
    """The GO DAG.

    Basically a dictionary between identifiers (i.e., strings like
    ``"GO:0006468"``) and ``GOTerm``'s.

    .. todo::

        Double check the depth computation, it seems to be a bit off.

    :param path: path to the OBO file.
    :param keep_obsolete: whether to keep obsolete terms (default: ``True``).
    """
    def __init__(self, path = None, keep_obsolete = True):
        self._terms = {}
        self._path = path
        if not path:
            return
        for term in _GOReader(path, keep_obsolete = keep_obsolete):
            assert not term.id in self._terms
            self._terms[term.id] = term
            for alt_id in term.alt_ids:
                assert not alt_id in self._terms
                self._terms[alt_id] = term
        self._populate_terms()

    def _set_term_depth(self, term):
        if term.level < 0:
            if not term.parents:
                term.level = 0
            else:
                term.level = min(self._set_term_depth(self._terms[parent])
                                 for parent, relation in term.parents
                                 if relation == "is_a") + 1
        return term.level

    def _populate_terms(self):
        for term in self._terms.itervalues():
            for parent, relation in term.parents:
                self._terms[parent].children.append((term.id, relation))
            if term.level < 0:
                self._set_term_depth(term)

    def __repr__(self):
        return "GODag('{}')".format(self._path)

    def __str__(self):
        return "\n".join("{}: {}".format(id_, term)
                         for id_, term in self._terms.items())

    def __getitem__(self, term_id):
        """Returns the ``GOTerm`` associated to a term ``id``."""
        if not term_id in self._terms:
            return
        return self._terms[term_id]

    def paths_to_root(self, term_id):
        """Returns all possible paths to the root node.

        Each path includes the term given. The order of the path is
        top -> bottom, i.e. it starts with the root and ends with the
        given term (inclusively).

        :param term_id: id of the GO term (e.g., ``'GO:0003682'``).
        :returns: a list of GO terms.
        """
        if term_id not in self._terms:
            return

        def _paths_to_root(term):
            if term.level == 0:
                return [[term]]
            paths = []
            for parent_id, relation in term.parents:
                paths_above = _paths_to_root(self._terms[parent_id])
                for path_above in paths_above:
                    paths.append([ term ] + path_above)
            return paths

        return _paths_to_root(self._terms[term_id])

    def get_valid_term_ids(self, include_bins=False):
        for id_ in self._terms.iterkeys():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_

    def get_valid_terms(self, include_bins=False):
        for id_, term in self._terms.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield term

    def get_valid_ids_terms(self, include_bins=False):
        for id_, term in self._terms.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_, term

    def get_interaspect_links(self):
        for term in self._terms.itervalues():
            for parent_id, relation in term.parents:
                parent = self._terms[parent_id]
                if term.namespace != parent.namespace:
                    yield term, relation, parent

    def add_bin_term(self, parent, proteins = set()):
        term = GOTerm()
        term.id = "BIN" + parent.id[parent_term.id.index(":"):]
        term.name = "Bin term for " + parent.name
        term.namespace = parent.namespace
        term.parents = [(parent.id, "bin_for")]
        term.children = []
        term.is_obsolete = False
        term.proteins = proteins
        term.level = parent.level + 1
        return term

    def _fill(self, p_to_terms):
        """Fills the DAG with proteins."""
        for p, terms in p_to_terms.iteritems():
            for term in terms:
                assert term.id in self._terms, "unknown term ID '{}'".format(term.id)
                term.proteins.add(p)

    def _generate_terms_by_level(self):
        """Generates all terms according to their recorder level."""
        level_to_terms = defaultdict(set)
        for term in self._terms.itervalues():
            level_to_terms[term.level].add(term)
        levels = sorted(level_to_terms.keys())
        for level in levels:
            for term in level_to_terms[level]:
                yield term

    def preprocess(self, p_to_terms, aspects = None, max_depth = None,
                   min_annot = None, add_bins = True):
        """Processes the DAG by removing unwanted terms and adding bin terms.

        This method does the following:

        * Removes from terms in unwanted namespaces.
        * Removes from terms with too few annotations.
        * Removes from terms that are too deep.
        * Adds *bin* terms.

        The protein list and protein-to-function map are adjusted accordingly.

        :param p_to_terms: map from protein IDs to ``GOTerm``'s.
        :param aspects: collection of aspects to keep; valid aspects are ``"biological_process"``, ``"cellular_component"``, ``"molecular_function"``.
        :param min_proteins_per_term: minimum number of proteins to keep a term.
        :param max_level: maximum term depth.
        :returns: the trimmed ``p_to_terms`` map, updates ``self``.
        """
        if aspects is None:
            aspects = ["biological_process", "cellular_component", "molecular_function"]

        self._fill(p_to_terms)

        for term in self._terms.itervalues():

            # Sanity check: check that a term has at most as many annotations
            # as the number of annotations in all its parents (note that a term
            # may have more than one parent, hence the sum).
            num_parent_annot = 0
            for parent_id, relation in term.parents:
                parent = self._terms[parent_id]
                if relation == "is_a":
                    num_parent_annot += len(parent.proteins)
            assert len(term.parents) == 0 or num_parent_annot >= len(term.proteins), \
                "failed sanity check: {} -> parents {}".format(term, term.parents)

            # Sanity check: check that a term has at least as many annotations
            # as each of its children
            for child_id, relation in term.children:
                child = self._terms[child_id]
                if relation == "is_a":
                    assert len(term.proteins) >= len(child.proteins), \
                        "failed sanity check: {} -> children {}".format(term, term.children)

        # Do a per-level traversal of the dag, and mark the terms that satisfy
        # all constraints. There may be cases where a term T has two child
        # terms C1 and C2 such that C2 is a child of C1 (i.e. C2 is a child of
        # both T and C1). DAGs be damned!
        processed, to_keep = set(), set()
        for term in self._generate_terms_by_level():
            if term in processed:
                continue
            processed.add(term)

            if not max_depth is None and term.level > max_depth:
                print "discarding '{}', too deep ({} > {})".format(term, term.level, max_depth)
                # Since this is a per-level traversal, we can break here
                break
            if not term.namespace in aspects:
                print "discarding '{}', not in namespace".format(term)
                continue
            if not min_annot is None and len(term.proteins) < min_annot:
                print "discarding '{}', too few annotations ({} < {})".format(term, len(term.proteins), min_annot)
                continue
            to_keep.add(term)

        # Remove all terms that are not in to_keep. For each term, if the
        # number of annotations in children terms that are marked for removal
        # is above min_annot, then add a bin node.
        # XXX how to deal with bin terms that are super-classes of terms that
        # have a common children with a term that is not marked for removal?
        raise NotImplementedError

        filtered_p_to_terms = {p: term for p, term in p_to_terms.iteritems()
                               if term in self._terms}
        return filtered_p_to_terms
