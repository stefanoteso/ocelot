.. ocelot documentation master file, created by
   sphinx-quickstart on Thu Feb  5 17:46:48 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Ocelot
======

Introduction
------------

Ocelot is a framework for building predictors of protein function and
interaction. It provides some python infrastructure to:

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

Ocelot is released under the XXX License. Components being used by Ocelot may
have widely different licensing constraints.

Requirements
------------

* Virtuoso
* rdflib
* SPARQLWrapper
* BioPython (optional)
* mmLib (optional)
* SBR
* SHOGUN (optional)

Usage
-----

* Make sure all the requirements (packages and data) are there.
* Run ``main.py make-rdf -s $PATH_TO_DB_DIR -d $PATH_TO_RDF_DIR``
* Run ``main.py upload-rdf -s $PATH_TO_RDF_DIR``
* Run ``main.py run-experiment -s $PATH_TO_DB_DIR -d $PATH_TO_RESULTS``

Authors
-------

- Stefano Teso (``name.surname _AT_ gmail.com``)
- Luca Masera


Contents
========

.. toctree::
   :maxdepth: 2

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
.. automodule:: ocelot.features
    :members:
.. automodule:: ocelot.kernels
    :members:
.. automodule:: ocelot.experiments
    :members:


References
==========

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

