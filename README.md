# entanglr
#### Bidirectional conversion between Seurat and AnnData (H5AD) formats with preservation of dimensionality reductions, neighbor graphs, and metadata. Built using  [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html).

## Overview

entanglr converts Seurat single-cell objects to AnnData H5AD format and vice versa. Unlike existing conversion tools, it preserves computed dimensional reductions (PCA, UMAP, t-SNE), neighbor graphs, PCA loadings, and variable feature annotations.

## Installation
entanglr can be installed from GitHub
```
devtools::install_github("neuroconnell/entanglr")
```

## Load the required packages
```
library(entanglr)
library(Seurat)
```

## Create Seurat object
```
seurat_obj <- readRDS("data.rds")
```

## Convert to H5AD
```
writeH5AD(
  seurat_object = seurat_obj,
  file = "output.h5ad",
  assay = "RNA",
  X_name = "data",
  raw = TRUE,
  reductions = TRUE,
  graphs = TRUE,
  var_features = TRUE,
  verbose = TRUE
  )
```
## Read H5AD and create Seurat object
```
seurat_obj <- readH5AD(
  file = "input.h5ad",
  assay = "RNA",
  verbose = TRUE
  )
```
## Dimensional reductions and graphs are automatically restored
```
names(seurat_obj@reductions)
names(seurat_obj@graphs)
```
