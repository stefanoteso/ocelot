# Ocelot

## What is it?

Ocelot is a collection of Python classes and methods for processing biological
data and building predictors.

Currently it focuses on protein function and interactions; interfaces to
more kinds of data will be added in the future.

The Ocelot infrastructure provides support for:

- Converting data from a handful of databases into RDF format.
- Building prediction datasets from the RDF data by querying a local SPARQL endpoint.
- Computing biologically meaningful features and kernels.
- Running complex prediction experiments.

## License

This code is licensed under the the 3-Clause BSD License.

## Requirements

The following software must be installed:

- [OpenLink Virtuoso open source](https://github.com/openlink/virtuoso-opensource), tested with version 7.2; the open source version is quite enough.

The following Python packages must be installed:

- [NumPy](http://www.numpy.org/), tested with version 1.8.2
- [SciPy](http://www.scipy.org/), tested with version 0.14.1
- [Pandas](http://pandas.pydata.org/), tested with version 0.17.1
- [scikit-learn](http://scikit-learn.org/stable/), tested with version 0.15.2
- [Matplotlib](http://matplotlib.sourceforge.net/), tested with version 1.4.2
- [BioPython](http://biopython.org), tested with version 1.64
- [rdflib](https://rdflib.readthedocs.org/en/latest/), tested with version 4.1.2
- [SPARQLWrapper](https://pypi.python.org/pypi/SPARQLWrapper), tested with version 1.6.0
- [pytest](http://pytest.org/latest/)


## How to use

Several settings are controlled by the ``ocelot.ini`` file.

- The first step is to convert the target biological datasets to RDF form.
  Make sure that the data is where the code can find it

```bash
$ path ./main.py make-rdf -s $path-to-root-data-directory -d $path-to-destination-directory
```

Please see the ``main.py::_make_rdf`` function for the structure that the
source directory should have.  The resulting RDF file will be stored in the
destination directory.

- Next you'll want to start up your ``virtuoso`` server and load the RDF data.
  Make sure that the code knows where to find the ``virtuoso`` executable.
This is controlled by ``ocelot.ini`` but can be overridden from the command
line.  Type ``./main.py --help`` to see the list of options.

```bash
$ python ./main.py start-virtuoso
$ python ./main.py upload-rdf
```

This takes care of uploading the RDF file produced above into a Virtuoso
knowledge graph.  The destination is controlled, again, by ``ocelot.ini``.

- Once the knowledge graph is in Virtuoso and can be efficiently queried, you
  can actually produce the annotations, rules, and kernels required by
SBR---already split into folds.  To produce the exact data used in the paper,
run:

```
$ python ./main.py run-experiment sgd
```

The resulting files will be stored in the destination directory specified
by the ``-d`` option.  These files can be used directly with SBR to perform
cross-validation.

## Authors

- Stefano Teso
- Luca Masera
