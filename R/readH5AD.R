#' Read H5AD File and Convert to Seurat with Full Restoration
#'
#' Reads an AnnData H5AD file and converts to Seurat, restoring dimensional
#' reductions, neighbor graphs, variable features, and scaled data if present.
#'
#' @param file Path to H5AD file
#' @param assay Name for the assay in Seurat object (default: "RNA")
#' @param verbose Print progress (default: TRUE)
#' @return A Seurat object
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- readH5AD("data.h5ad")
#' }
readH5AD <- function(file, assay = "RNA", verbose = TRUE) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  if (verbose) message("Reading H5AD file...")

  .h5ad_to_list <- function(file) {
    anndata <- reticulate::import("anndata")

    adata <- anndata$read_h5ad(file)

    result <- list(
      X = as.matrix(adata$X),
      obs = as.data.frame(adata$obs),
      var = as.data.frame(adata$var),
      raw = NULL,
      obsm = list(),
      varm = list(),
      obsp = list(),
      layers = list()
    )

    # Get raw counts if available
    if (!is.null(adata$raw)) {
      result$raw <- as.matrix(adata$raw$X)
    }

    # Get dimensional reductions from obsm
    tryCatch({
      pca_data <- adata$obsm[["X_pca"]]
      if (!is.null(pca_data)) {
        result$obsm[["X_pca"]] <- as.matrix(pca_data)
      }
    }, error = function(e) {})

    tryCatch({
      umap_data <- adata$obsm[["X_umap"]]
      if (!is.null(umap_data)) {
        result$obsm[["X_umap"]] <- as.matrix(umap_data)
      }
    }, error = function(e) {})

    tryCatch({
      tsne_data <- adata$obsm[["X_tsne"]]
      if (!is.null(tsne_data)) {
        result$obsm[["X_tsne"]] <- as.matrix(tsne_data)
      }
    }, error = function(e) {})

    # Get loadings from varm (e.g., PCA loadings)
    tryCatch({
      pcs <- adata$varm[["PCs"]]
      if (!is.null(pcs)) {
        result$varm[["PCs"]] <- as.matrix(pcs)
      }
    }, error = function(e) {})

    # Get graphs from obsp
    tryCatch({
      conn <- adata$obsp[["connectivities"]]
      if (!is.null(conn)) {
        result$obsp[["connectivities"]] <- as.matrix(conn)
      }
    }, error = function(e) {})

    tryCatch({
      dist <- adata$obsp[["distances"]]
      if (!is.null(dist)) {
        result$obsp[["distances"]] <- as.matrix(dist)
      }
    }, error = function(e) {})

    # Get additional layers
    tryCatch({
      scaled <- adata$layers[["scaled"]]
      if (!is.null(scaled)) {
        result$layers[["scaled"]] <- as.matrix(scaled)
      }
    }, error = function(e) {})

    return(result)
  }

  adata_list <- basilisk::basiliskRun(
    env = entanglr_env,
    fun = .h5ad_to_list,
    file = normalizePath(file)
  )

  if (verbose) message("Converting to Seurat...")

  # Create count matrix
  counts <- if (!is.null(adata_list$raw)) {
    if (verbose) message("Using raw counts from .raw slot")
    Matrix::Matrix(t(adata_list$raw), sparse = TRUE)
  } else {
    if (verbose) message("No raw counts found, using X as counts")
    Matrix::Matrix(t(adata_list$X), sparse = TRUE)
  }

  if (nrow(adata_list$var) > 0) {
    rownames(counts) <- rownames(adata_list$var)
  }
  if (nrow(adata_list$obs) > 0) {
    colnames(counts) <- rownames(adata_list$obs)
  }

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay,
    meta.data = adata_list$obs
  )

  # Add normalized data if different from raw
  if (!is.null(adata_list$raw)) {
    if (verbose) message("Adding normalized data from X")
    norm_data <- Matrix::Matrix(t(adata_list$X), sparse = TRUE)

    if (nrow(adata_list$var) > 0) {
      rownames(norm_data) <- rownames(adata_list$var)
    }
    if (nrow(adata_list$obs) > 0) {
      colnames(norm_data) <- rownames(adata_list$obs)
    }

    seurat_obj <- tryCatch({
      Seurat::SetAssayData(
        seurat_obj,
        layer = "data",
        new.data = norm_data
      )
    }, error = function(e) {
      Seurat::SetAssayData(
        seurat_obj,
        slot = "data",
        new.data = norm_data
      )
    })
  }

  # Restore dimensional reductions
  if (length(adata_list$obsm) > 0) {
    if (verbose) message("Restoring dimensional reductions...")

    for (obsm_key in names(adata_list$obsm)) {
      tryCatch({
        embeddings <- adata_list$obsm[[obsm_key]]

        # Map back to Seurat names
        reduction_name <- gsub("^X_", "", obsm_key)

        # Set rownames to match cells
        if (ncol(seurat_obj) == nrow(embeddings)) {
          rownames(embeddings) <- colnames(seurat_obj)
        }

        # Create key matching Seurat conventions
        key <- switch(reduction_name,
                      "pca" = "PC_",
                      "umap" = "UMAP_",
                      "tsne" = "tSNE_",
                      paste0(toupper(substr(reduction_name, 1, 1)),
                             substr(reduction_name, 2, nchar(reduction_name)), "_"))

        # Set column names to match original format
        colnames(embeddings) <- paste0(key, 1:ncol(embeddings))

        # Get or create loadings if available
        loadings <- NULL
        if (reduction_name == "pca" && "PCs" %in% names(adata_list$varm)) {
          loadings <- adata_list$varm[["PCs"]]
          # Ensure rownames match genes
          if (nrow(seurat_obj) >= nrow(loadings)) {
            matching_genes <- intersect(rownames(seurat_obj), rownames(loadings))
            if (length(matching_genes) > 0) {
              loadings <- loadings[matching_genes, , drop = FALSE]
              # Set column names for loadings
              colnames(loadings) <- paste0(key, 1:ncol(loadings))
            }
          }
        }

        # Create DimReduc object
        if (!is.null(loadings) && nrow(loadings) > 0) {
          seurat_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
            embeddings = embeddings,
            loadings = loadings,
            key = key,
            assay = assay
          )
        } else {
          seurat_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
            embeddings = embeddings,
            key = key,
            assay = assay
          )
        }

        if (verbose) message(sprintf("  - Restored %s", reduction_name))

      }, error = function(e) {
        if (verbose) warning(sprintf("Could not restore %s: %s", obsm_key, e$message))
      })
    }
  }

  # Restore neighbor graphs
  if (length(adata_list$obsp) > 0) {
    if (verbose) message("Restoring neighbor graphs...")

    for (obsp_key in names(adata_list$obsp)) {
      tryCatch({
        graph_matrix <- adata_list$obsp[[obsp_key]]

        # Convert to sparse matrix and transpose back
        graph <- methods::as(Matrix::t(graph_matrix), "dgCMatrix")

        # Set row and column names
        if (ncol(seurat_obj) == nrow(graph)) {
          rownames(graph) <- colnames(seurat_obj)
          colnames(graph) <- colnames(seurat_obj)
        }

        # Map back to Seurat naming
        graph_name <- switch(obsp_key,
                             "connectivities" = paste0(assay, "_snn"),
                             "distances" = paste0(assay, "_nn"),
                             obsp_key)

        seurat_obj@graphs[[graph_name]] <- graph

        if (verbose) message(sprintf("  - Restored %s graph", obsp_key))

      }, error = function(e) {
        if (verbose) warning(sprintf("Could not restore graph %s: %s", obsp_key, e$message))
      })
    }
  }

  # Restore variable features
  if ("highly_variable" %in% colnames(adata_list$var)) {
    tryCatch({
      var_genes <- rownames(adata_list$var)[adata_list$var$highly_variable]
      if (length(var_genes) > 0) {
        Seurat::VariableFeatures(seurat_obj, assay = assay) <- var_genes
        if (verbose) message(sprintf("Restored %d variable features", length(var_genes)))
      }
    }, error = function(e) {
      if (verbose) warning(sprintf("Could not restore variable features: %s", e$message))
    })
  }

  # Restore scaled data
  if ("scaled" %in% names(adata_list$layers)) {
    tryCatch({
      scaled_matrix <- adata_list$layers[["scaled"]]
      scaled <- methods::as(Matrix::t(scaled_matrix), "dgCMatrix")

      # Set dimnames
      if (nrow(adata_list$var) > 0) {
        rownames(scaled) <- rownames(adata_list$var)
      }
      if (nrow(adata_list$obs) > 0) {
        colnames(scaled) <- rownames(adata_list$obs)
      }

      seurat_obj <- tryCatch({
        Seurat::SetAssayData(
          seurat_obj,
          layer = "scale.data",
          new.data = scaled
        )
      }, error = function(e) {
        Seurat::SetAssayData(
          seurat_obj,
          slot = "scale.data",
          new.data = scaled
        )
      })

      if (verbose) message("Restored scaled data")

    }, error = function(e) {
      if (verbose) message("Could not restore scaled data")
    })
  }

  # Restore active identity - check common column names
  if ("tree.ident" %in% colnames(seurat_obj@meta.data)) {
    if (verbose) message("Setting active identity to 'tree.ident'")
    Idents(seurat_obj) <- seurat_obj$tree.ident
  } else if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    if (verbose) message("Setting active identity to 'seurat_clusters'")
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
  } else if ("ident" %in% colnames(seurat_obj@meta.data)) {
    if (verbose) message("Setting active identity to 'ident'")
    Idents(seurat_obj) <- seurat_obj$ident
  }

  if (verbose) message("Successfully converted to Seurat object")
  return(seurat_obj)
}
