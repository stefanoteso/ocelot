# Ocelot

![Ocelot Icon](https://upload.wikimedia.org/wikipedia/commons/2/27/Salvador_Dali_NYWTS.jpg)

## What is it?

Ocelot is a collection of Python utility classes and methods for processing
biological data and building predictors.

Currently it focuses on protein function and interactions; interfaces to
more kinds of data will be added as they are required.

The Ocelot infrastructure provides support for:

-   Converting data from a handful of databases into RDF format.
-   Building prediction datasets from the RDF data by querying a local SPARQL
    endpoint.
-   Computing biologically meaningful features and kernels.
-   Running complex prediction experiments.

## License

WRITEME

## Requirements

-   [OpenLink Virtuoso open source](https://github.com/openlink/virtuoso-opensource), tested with version 7.2
-   [rdflib](https://rdflib.readthedocs.org/en/latest/), tested with version 4.1.2
-   [SPARQLWrapper](https://pypi.python.org/pypi/SPARQLWrapper), tested with version 1.6.0
-   [NumPy](http://www.numpy.org/), tested with version 1.8.2
-   [SciPy](http://www.scipy.org/), tested with version 0.14.1
-   [scikit-learn](http://scikit-learn.org/stable/), tested with version 0.15.2
-   [SHOGUN](http://www.shogun-toolbox.org/) through the modular Python interface, tested with version 3.2.1
-   [Semantic Based Regularization](https://sites.google.com/site/semanticbasedregularization/home/software), tested with version 1.1.1
-   [Matplotlib](http://matplotlib.sourceforge.net/), tested with version 1.4.2
-   [BioPython](http://biopython.org), tested with version 1.64
-   [mmLib](http://pymmlib.sourceforge.net/)
-   [pytest](http://pytest.org/latest/)

## How to use

Make sure that you have ``virtuoso`` installed.

-   The first step is to convert the target biological datasets to RDF form.

    WRITEME

-   Next you'll want to start up your ``virtuoso`` server and load the RDF data.

    WRITEME

-   Then you'll want to run your experiment.

    WRITEME

## Authors

- Stefano Teso
- Luca Masera
