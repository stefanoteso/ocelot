# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from ocelot.kernels import Kernel

def read_pcl(path, id_column="NAME"):
    """Reads a PCL microarray file [1]_.

    The file is simply a TSV where the first three columns are fixed,
    followed by a variable number of columns (one per condition).

    It currently ignores GWEIGHT and EWEIGHT.

    Parameters
    ----------
    path : str or any object with a read() method.
        Path to the PCL file.
    id_column : str, optional. (defaults to ``"NAME"``)
        Name of the column to use for gene IDs.

    Returns
    -------
    data : pandas.DataFrame
        Map from gene ID to expression levels (list of floats).

    References
    ----------
    .. [1] https://puma.princeton.edu/help/formats.shtml#pcl
    """
    df = pd.read_csv(path, sep="\t", header=0)[1:]
    if df.isnull().values.any():
        raise ValueError("PCL file has missing entries")
    df = df.set_index(id_column)
    df = df.drop(df.columns[0], axis=1)
    df = df.drop(df.columns[0], axis=1)
    return df

class GeneExprKernel(Kernel):
    """A kernel on microarray data.

    Parameters
    ----------
    gene_to_i : dict
        Map from gene names to indices in the Gram matrix. There should be
        a key for each target gene.
    microarrays : list
        List of microarray experiments, as returned by e.g. ``read_pcl()``.
    """
    def __init__(self, gene_to_i, microarrays, *args, **kwargs):
        self._microarrays = microarrays
        super(GeneExprKernel, self).__init__(gene_to_i, *args, **kwargs)

    def _compute_all(self):
        gene_to_i = self._entities
        target_genes, _ = zip(*sorted(gene_to_i.items(),
                                      key=lambda gene_i: gene_i[-1]))

        matrix = np.zeros((len(self), len(self)), dtype=self.dtype)
        for microarray in self._microarrays:
            # Sort the microarray so that the first N rows contain the
            # expression levels of target genes (in the correct order), and
            # the remaining rows contain non-target genes (in any order).
            target_genes_in_microarray = [gene for gene in target_genes
                                          if gene in microarray.index]
            other_genes_in_microarray = [gene for gene in microarray.index
                                         if gene not in target_genes]
            sorted_genes = target_genes_in_microarray + other_genes_in_microarray
            assert len(sorted_genes) == microarray.shape[0]
            microarray.reindex(sorted_genes)

            # Compute the correlation, fix cases with zero stddev
            corr = np.corrcoef(microarray.values)
            corr[np.isnan(corr).nonzero()] = 0.0

            # Add to the full matrix
            target_indices = [gene_to_i[gene] for gene
                              in target_genes_in_microarray]
            matrix[np.ix_(target_indices, target_indices)] += \
                corr[:len(target_genes_in_microarray),
                     :len(target_genes_in_microarray)]
        return matrix
