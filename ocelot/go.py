# -*- coding: utf-8 -*-

import sys, re

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
        """ID of the term, e.g. ``"GO:GO:0048308'"``"""
        self.alt_ids = []
        """List of alternative IDs."""
        self.name = ""
        """Short description of the term."""
        self.namespace = ""
        """Namespace of the term, one of ``["BP", "CC", "MF"]``."""
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

#    def get_all_parent_edges(self, dag):
#        all_parent_edges = set()
#        for p, r in self.parents:
#            all_parent_edges.add((self.id, dag[p].id, r))
#            all_parent_edges |= dag[p].get_all_parent_edges(dag)
#        return all_parent_edges

#    def get_all_child_edges(self, dag):
#        all_child_edges = set()
#        for p, r in self.children:
#            all_child_edges.add((dag[p].id, self.id, r))
#            all_child_edges |= dag[p].get_all_child_edges(dag)
#        return all_child_edges

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

    def add_bin_term(parent, proteins = set()):
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

    def get_valid_term_ids(self, include_bins=False):
        for id_ in self._terms.iterkeys():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_

    def get_valid_terms(self, include_bins=False):
        for id_, term in self.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield term

    def get_valid_id(self, include_bins=False):
        for id_, term in self.iteritems():
            if id_.startswith("GO:") or \
               (include_bins and id_.startswith("BN:")):
                yield id_, term

#    def prune(self, ids_to_keep, 
#        from copy import deepcopy
#
#        dag = GODag()
#        for id_ in id_to_keep:
#            dag._terms[id_] = deepcopy(self._terms[id_])
#
#        for id_, term in dag._terms.items():
#            term.proteins = set()
#
#        for key, t in new_dag.items():
#            t.proteins = set()
#            new_children = []
#            new_parents = []
#            proteins_in_the_bin = set()
#            for c, r in t.children:
#                if c in new_dag:
#                    new_children += [(c, r)]
#                else:
#                    if r == 'is_a' and c in proteins_by_go_id:
#                        proteins_in_the_bin.update(proteins_by_go_id[c])
#            if new_children:
#                if proteins_in_the_bin:
#                    bin_children = create_bin_node(t, proteins_in_the_bin if not valid_proteins else proteins_in_the_bin.intersection(valid_proteins))
#                    new_children += [(bin_children.id, 'bin')]
#                    new_dag[bin_children.id] = bin_children
#            else:
#                t.proteins = proteins_by_go_id[key] if not valid_proteins else proteins_by_go_id[key].intersection(valid_proteins)
#
#            t.children = new_children
#
#            for p, r in t.parents:
#                if p in new_dag:
#                    new_parents += [(p, r)]
#            t.parents = new_parents
#
#        for key, t in new_dag.iteritems():
#            if not t.children:
#                for p in t.get_all_parents(new_dag, rel=['is_a', 'bin', 'part_of']):
#                    new_dag[p].proteins.update(t.proteins)
#
#        return new_dag
