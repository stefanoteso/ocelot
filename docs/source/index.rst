.. ocelot documentation master file, created by
   sphinx-quickstart on Thu Feb  5 17:46:48 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Ocelot
======

Introduction
------------

Ocelot is a framework for applying Statistical-Relational Learning [SRL]_,
[Getoor07]_, [DeRaedt08]_ methods to the prediction of protein features and
interactions from either sequence or structure.

Current experiments are all based on Semantic Based Regularization
[Diligenti12]_.  Prior work on the same subject can be found in [Sacca14]_.

Requirements
------------

Depending on what you want to do with ocelot, the requirements change.

#. You will need a working installation of Virtuoso [Virtuoso]_.

Experiment
----------

Steps to reproduce the experiment:

#. Make sure all the requirements (packages and data) are there.
#. Run ``main.py make-rdf -s $PATH_TO_DATABASE_DIR -d $PATH_TO_RDF_DIR``
#. Run ``main.py upload-rdf -s $PATH_TO_RDF_DIR``
#. Run ``main.py run-experiment ...``


Contents
========

.. toctree::
   :maxdepth: 2

.. automodule:: ocelot.ontology
    :members:
.. automodule:: ocelot.converters.biogrid
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

.. [SRL] https://en.wikipedia.org/wiki/Statistical_relational_learning

.. [Getoor07] Getoor and Taskar, `Introduction to Statistical Relational Learning`, 2007, The MIT Press

.. [DeRaedt08] De Raedt and Kersting, `Probabilistic inductive logic programming', 2008, Springer

.. [Diligenti12] Diligenti et al., `Bridging Logic and Kernel Machines`, 2012, Machine Learning

.. [Sacca14] Sacca' et al., `Improved multi-level protein-protein interaction prediction with Semantic-based Regularization`, 2014, BMC Bioinformatics

.. [Virtuoso] http://virtuoso.openlinksw.com/dataspace/doc/dav/wiki/Main/


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

