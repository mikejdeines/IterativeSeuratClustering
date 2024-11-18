Performs Leiden clustering (Traag, V. A. et al. Sci Rep 9, 5233 (2019)) iteratively on Seurat objects to identify all clusters with a significant number of differentially-expressed genes.
The gene scoring metric was adapted from Tasic, B. et al. Nat Neurosci 19, 335–346 (2016), with minor modifications for differential expression testing and simplified cutoffs.

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
Due to some quirks in the SplitObject function, some of the downstream analysis needs to be modified:

RunUMAP should be performed prior to clustering.
If the data has been SCTransformed, FindMarkers or FindAllMarkers should be run with the parameter "recorrect_umi = FALSE".
