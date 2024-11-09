Performs Leiden clustering (Traag, V. A. et al. Sci Rep 9, 5233 (2019)) iteratively on Seurat objects to identify all clusters with a significant number of differentially-expressed genes.

Installation:

```
devtools::install_github("mikejdeines/IterativeSeuratClustering")
```

This package requires the use of igraph's implementation of the Leiden algorithm.
To use this implementation, create a conda environment using reticulate and install numpy 1.24.4, igraph, and leidenalg.

```
library(reticulate)
conda_create(envname = "scrnaseq")
use_condaenv("scrnaseq")
py_install("numpy==1.24.4")
py_install("igraph")
py_install("leidenalg")
```
