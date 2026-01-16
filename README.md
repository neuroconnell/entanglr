# entanglr
#### R package for bidirectional conversion between [Seurat](https://github.com/satijalab/seurat) and [AnnData](https://github.com/scverse/anndata) (H5AD) data formats while preserving dimensionality reductions, neighbor graphs, and metadata. Python dependencies are managed automatically via  [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html).

## Overview

**entanglr** directly converts Seurat objects to AnnData H5AD format and vice versa. It preserves computed dimensional reductions (PCA, UMAP, t-SNE), neighbor graphs, PCA loadings, and variable feature annotations.

### What additional data does entanglr maintain?

**Dimensional reductions:**
- Cell embeddings (PCA, UMAP, tSNE coordinates) with zero numerical loss
- PCA feature loadings (zero-padded to all genes)
- Reduction keys (PC_, UMAP_, tSNE_)
- **Support for multiple reductions** (e.g., `pca`, `integrated.pca`, `umap`, `integrated.umap`)

**Neighbor graphs:**
- Graph structure and edge weights with zero numerical loss
- Multiple graphs (SNN, KNN) from different assays
- **Graph names preserved** (e.g., `integrated_snn` remains `integrated_snn`)

**Variable features:**
- Highly variable genes
- **Gene order preserved** via ranking system
- Automatically extracted from PCA loadings if not explicitly stored

## Key Features

**Dynamic metadata reading** - Not limited to hardcoded reduction/graph names  

**Multi-assay support** - Preserves both RNA and integrated assay reductions  

**Zero data loss** - Values are preserved with no loss

-------

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

## Converting Seurat object to H5AD
Note: First use triggers automatic Python environment setup (~ 5 minutes). Writing of all subsequent H5AD files is immediate.

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
## Reading H5AD file to create a Seurat object
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










