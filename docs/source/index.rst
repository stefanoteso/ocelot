.. ocelot documentation master file, created by
   sphinx-quickstart on Thu Feb  5 17:46:48 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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

* `OpenLink Virtuoso open-source <https://github.com/openlink/virtuoso-opensource>`_
* `rdflib <https://rdflib.readthedocs.org/en/latest/>`_
* `SPARQLWrapper <https://pypi.python.org/pypi/SPARQLWrapper>`_
* Optionally `BioPython <http://biopython.org>`_
* Optionally `mmLib <http://pymmlib.sourceforge.net/>`_
* `Semantic Based Regularization <https://sites.google.com/site/semanticbasedregularization/home/software>`_
* `SHOGUN <http://www.shogun-toolbox.org/>`_
* `scikit-learn <http://scikit-learn.org/stable/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_

Usage
-----

First of all, make sure that all the `Requirements`_ are in place.

Assuming that you want to run the ``yip09`` experiment, you'll first need to
build the RDF dataset starting from the (various) database dumps. From the
command line, run::

    $ ./main.py make-rdf -s $PATH_TO_DATABASES -d rdf-data -t sgd,yip09

This will take the ``SGD`` and ``yip09`` dumps, turn them into RDF triples,
and place the result in the ``rdf-data`` directory.

Now you should load the RDF data into your Virtuoso instance of choice, using
your method of choice. Ocelot provides a shortcut for loading the data to a
*local* Virtuoso instance via the ``isql`` command. This can be done as
follows::

    $ ./main.py upload-rdf -s rdf-data -g "http://ocelot-yip09-graph"

where ``http://ocelot-yip09-graph`` is the RDF graph to be used.

Now you can run the experiment by typing the following::

    $ ./main.py run-experiment -s $PATH_TO_DATABASES -e "http://localhost:8890/sparql" -g "http://ocelot-yip09-graph"

where ``http://localhost:8890/sparql`` is the URI of the SPARQL endpoint of
your Virtuoso instance. The path to the database dumps should be provided here
as well, since some data (e.g. the PSSM files and the microarray files) is
*not* actually converted to RDF (and it would be pointless to do so).

Authors
-------

- Stefano Teso (``name.surname _AT_ gmail.com``)

References
----------

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

**Yip et al. experiment**

.. [Yip09] Yip, Kim, McDermott, Gerstein, *Multi-level learning: improving the
    prediction of protein, domain and residue interactions by allowing
    information flow between levels*, BMC Bioinformatics, 2009.

**Yeast Datasets**

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

**Kernels**

.. [Leslie02a] Leslie et al., *The spectrum kernel: A string kernel for SVM
    protein classification*, 2002.

.. [Leslie02b] Leslie et al., *Mismatch String Kernels for SVM Protein
    Classification*, 2002.

.. [Kuang04] Kuang et al., *Profile-based string kernels for remote homology
    detection and motif extraction*, 2004.

.. [Kondor02] Kondor and Lafferty, *Diffusion Kernels on Graphs and Other
    Discrete Input Spaces*, 2002.


Contents
========

.. toctree::
   :maxdepth: 2

**RDF**

.. automodule:: ocelot.ontology
    :members:

**Services**

.. automodule:: ocelot.services
    :members:

**Converters**

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

**Features and Kernels**

.. automodule:: ocelot.features
    :members:
.. automodule:: ocelot.kernels.base
    :members:
.. automodule:: ocelot.kernels.vector
    :members:
.. automodule:: ocelot.kernels.string
    :members:
.. automodule:: ocelot.kernels.graph
    :members:

**Experiments**

.. automodule:: ocelot.experiments
    :members:




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

