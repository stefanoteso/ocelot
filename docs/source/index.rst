Ocelot
======

Ocelot is a collection of scripts for building predictors of protein function
and interaction. It provides some Python infrastructure to:

#. Convert some databases into RDF format (though currently it does not
   make use of a formal TBox), serialized in turtle notation.
#. Build prediction datasets from the RDF data by querying a local SPARQL
   endpoint.
#. Computing biologically significant features and kernels, exploiting
   information from protein/domain/residue sequence, homology and structure,
   with additional data sources such as gene expression data, protein
   complexes, and domain annotations.
#. Run protein function/interactions experiments.

Ocelot is being used for studying the applicability of Statistical-Relational
Learning [SRL]_, [Getoor07]_, [DeRaedt08]_ methods to predictive proteomics.
The current experiments are all based on Semantic Based Regularization
[Diligenti12]_; some prior work on the subject can be found in [Sacca14]_.

License
-------

WRITEME

Requirements
------------

For dealing with RDF:

* `OpenLink Virtuoso open-source <https://github.com/openlink/virtuoso-opensource>`_, tested with version ``7.2.1``
* `rdflib <https://rdflib.readthedocs.org/en/latest/>`_, tested with version ``4.1.2``
* `SPARQLWrapper <https://pypi.python.org/pypi/SPARQLWrapper>`_, tested with version ``1.6.0``

For dealing with other data:

* `NumPy <http://www.numpy.org/>`_, tested with version ``1.8.2``
* `SciPy <http://www.scipy.org/>`_, tested with version ``0.14.1``
* Optionally `BioPython <http://biopython.org>`_, tested with version ``1.64``
* Optionally `mmLib <http://pymmlib.sourceforge.net/>`_

For the learning algorithms:

* `Semantic Based Regularization <https://sites.google.com/site/semanticbasedregularization/home/software>`_, tested with version ``1.1.1``
* `SHOGUN <http://www.shogun-toolbox.org/>`_ (modular python interface), tested with version ``3.2.1``
* `scikit-learn <http://scikit-learn.org/stable/>`_, tested with version ``0.15.2``
* Optionally `Matplotlib <http://matplotlib.sourceforge.net/>`_, tested with version ``1.4.2``

To build the documentation:

* `sphinx <http://sphinx-doc.org/>`_, tested with version ``1.2.3``, as well as the `todo plugin <http://sphinx-doc.org/ext/todo.html>`_

To run the (very few) tests:

* `pytest <http://pytest.org/latest/>`_

Usage
-----

First of all, make sure that all the `Requirements`_ are in place.

You can glance at the Ocelot help text to get an overview of the
commands it accepts and the available options::

    $ ./main.py --help

Assuming that you want to run the ``yip09`` experiment, you'll first need to
build the RDF dataset starting from the (various) database dumps.
Unfortunately, this stel is not automated (yet), so you will need to manually
fetch the required source databases. After you collected all of them, run the
following::

    $ ./main.py make-rdf -s $PATH_TO_DATABASES -d rdf-data -t sgd,yip09

This will take the ``SGD`` and ``yip09`` dumps, turn them into RDF triples,
and place the result in the ``rdf-data`` directory.

Now you can (and should) load the RDF data into your Virtuoso instance of
choice, using your method of choice. Ocelot provides a shortcut for loading the
data to a *local* Virtuoso instance via the ``isql`` command. This can be done
as follows::

    $ ./main.py upload-rdf -s rdf-data -g "http://ocelot-yip09-graph"

where ``http://ocelot-yip09-graph`` is the RDF graph to be used. The ``-g``
option is *optional*; by default Ocelot will use the default graph specified
in the ``ocelot.ini`` configuration file.

Now you can run the experiment by typing the following::

    $ ./main.py run-experiment -t yip09 -s $PATH_TO_DATABASES -e "http://localhost:8890/sparql" -g "http://ocelot-yip09-graph"

where ``http://localhost:8890/sparql`` is the URI of the SPARQL endpoint of
your Virtuoso instance. The path to the database dumps should be provided here
as well, since some data (e.g. the PSSM files and the microarray files) is
*not* actually converted to RDF (and it would be pointless to do so).

In case you want to get rid of a graph stored in your local Virtuoso instance,
Ocelot provides a shortcut (implemented via the ``isql`` command)::

    $ ./main.py clear-graph -g "http://ocelot-yip09-graph

You can run the (admittedly very thin) test suite by running::

    $ py.test ocelot

from the main Ocelot directory.

Authors
-------

- Stefano Teso (``name.surname _AT_ gmail.com``)
- Luca Masera (WRITEME)


Contents
========

.. toctree::
   :maxdepth: 2

Converters
----------

.. automodule:: ocelot.ontology
    :members:
.. automodule:: ocelot.converters.ipfam
    :members:
.. automodule:: ocelot.converters.pdb
    :members:
.. automodule:: ocelot.converters.psimi
    :members:
.. automodule:: ocelot.converters.sgd
    :members:
.. automodule:: ocelot.converters.sifts
    :members:
.. automodule:: ocelot.converters.string
    :members:
.. automodule:: ocelot.converters.yip09
    :members:

Services
--------

.. automodule:: ocelot.services
    :members:
.. automodule:: ocelot.go
    :members:


Features
--------

.. automodule:: ocelot.features
    :members:

Kernels
-------

.. automodule:: ocelot.kernels
    :members:
.. automodule:: ocelot.kernels.vector
    :members:
.. automodule:: ocelot.kernels.string
    :members:
.. automodule:: ocelot.kernels.graph
    :members:
.. automodule:: ocelot.kernels.other
    :members:

Experiments
-----------

.. automodule:: experiments.yip09
    :members:
.. automodule:: experiments.yeast
    :members:


References
==========

Machine Learning
----------------

.. [SRL] `<https://en.wikipedia.org/wiki/Statistical_relational_learning>`_

.. [Getoor07] Getoor and Taskar, *Introduction to Statistical Relational
    Learning*, 2007, The MIT Press

.. [DeRaedt08] De Raedt and Kersting, *Probabilistic inductive logic
    programming*, 2008, Springer

.. [Diligenti12] Diligenti et al., *Bridging Logic and Kernel Machines*, 2012,
    Machine Learning

.. [Sacca14] Sacca' et al., *Improved multi-level protein-protein interaction
    prediction with Semantic-based Regularization*, 2014, BMC Bioinformatics

.. [Virtuoso] `<http://virtuoso.openlinksw.com/dataspace/doc/dav/wiki/Main/>`_

Features
--------

.. [WikipediaAAa] https://en.wikipedia.org/wiki/Amino_acid

.. [WikipediaAAb] https://en.wikipedia.org/wiki/Proteinogenic_amino_acid

.. [Swart04] Swart et al., *Polarizabilities of amino acid residues*, 2004.

.. [Wootton94] Wootton, *Sequences with unusual amino acid compositions*,
    1994.

.. [Liao05] Liao et al., *Protein sequence entropy is closely related to
    packing density and hydrophobicity*, 2005.

.. [Coletta10] Coletta et al., *Low-complexity regions within protein sequences
    have position-dependent roles*, 2010.

.. [Scholkopf99] Scholkopf et al., *Input Space Versus Feature Space in
    Kernel-Based Methods*, 1999.

Kernels
-------

.. [Leslie02a] Leslie et al., *The spectrum kernel: A string kernel for SVM
    protein classification*, 2002.

.. [Leslie02b] Leslie et al., *Mismatch String Kernels for SVM Protein
    Classification*, 2002.

.. [Kuang04] Kuang et al., *Profile-based string kernels for remote homology
    detection and motif extraction*, 2004.

.. [Kondor02] Kondor and Lafferty, *Diffusion Kernels on Graphs and Other
    Discrete Input Spaces*, 2002.

Yip09/Yeast Experiments
-----------------------

.. [Yip09] Yip, Kim, McDermott, Gerstein, *Multi-level learning: improving the
    prediction of protein, domain and residue interactions by allowing
    information flow between levels*, BMC Bioinformatics, 2009.

.. [Ito00] Ito et al. *Toward a Protein-Protein Interaction Map of the Budding
    Yeast: A Comprehensive System to Examine Two-Hybrid Interactions in All
    Possible Combinations between the Yeast Proteins*, PNAS, 2000.

.. [Uetz00] Uetz et al., *A Comprehensive Analysis of Protein-Protein
    Interactions in Saccharomyces cerevisiae*, Nature, 2000.

.. [Gavin06] Gavin et al., *Proteome Survey Reveals Modularity of the Yeast
    Cell Machinery*, Nature, 2006.

.. [Krogan06] Krogan et al., *Global Landscape of Protein Complexes in the
    Yeast Saccharomyces cerevisiae*, Nature, 2006

.. [Pu08] Pu et al., *Up-to-date catalogues of yeast protein complexes*, NAR
    2008

.. [Lee03] Lee and Sonnhammer, *Genomic gene clustering analysis of pathways in
    eukaryotes*, 2003.

.. [Yu10] Yu et al., *Simple sequence-based kernels do not predict
    protein–protein interactions*, 2010.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

