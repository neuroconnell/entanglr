#' Convert Seurat Object to H5AD Format with Full Preservation
#'
#' Converts a Seurat object to AnnData H5AD format, preserving dimensional
#' reductions, neighbor graphs, variable features, and scaled data.
#'
#' @param seurat_object A Seurat object
#' @param file Output file path (e.g., "data.h5ad")
#' @param assay Which assay to use (default: "RNA")
#' @param X_name Which data slot/layer to use (default: "data")
#' @param raw Include raw counts (default: TRUE)
#' @param reductions Include dimensional reductions like PCA, UMAP (default: TRUE)
#' @param graphs Include neighbor graphs (default: TRUE)
#' @param var_features Include variable features (default: TRUE)
#' @param scaled_data Include scaled data (default: FALSE)
#' @param verbose Print progress (default: TRUE)
#' @return Invisibly returns the file path
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' writeH5AD(seurat_obj, "output.h5ad", assay = "RNA")
#' }
writeH5AD <- function(seurat_object, file, assay = "RNA",
                      X_name = "data", raw = TRUE,
                      reductions = TRUE, graphs = TRUE,
                      var_features = TRUE, scaled_data = FALSE,
                      verbose = TRUE) {

  if (verbose) message("Using assay: ", assay)

  .seurat_to_h5ad <- function(seurat_object, file, assay, X_name, raw,
                              reductions, graphs, var_features, scaled_data, verbose) {

    anndata <- reticulate::import("anndata", convert = FALSE)
    scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
    pandas <- reticulate::import("pandas", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)

    if (verbose) message("Extracting Seurat data...")

    # Extract main data matrix
    X_data <- tryCatch({
      SeuratObject::GetAssayData(seurat_object, assay = assay, layer = X_name)
    }, error = function(e) {
      SeuratObject::GetAssayData(seurat_object, assay = assay, slot = X_name)
    })

    if (!inherits(X_data, "dgCMatrix")) {
      X_data <- methods::as(X_data, "dgCMatrix")
    }

    X_t <- Matrix::t(X_data)
    X_py <- reticulate::r_to_py(as.matrix(X_t))
    X_sparse <- scipy_sparse$csr_matrix(X_py)

    # Get metadata and preserve active identity
    obs <- seurat_object@meta.data

    # Save active identity if not already in metadata
    if (!"ident" %in% colnames(obs)) {
      obs$ident <- as.character(Idents(seurat_object))
    }

    var <- data.frame(row.names = rownames(X_data))

    # Mark variable features in var BEFORE creating dataframe
    if (var_features) {
      tryCatch({
        var_genes <- Seurat::VariableFeatures(seurat_object, assay = assay)

        # If no variable features but PCA exists, extract from PCA loadings
        if (length(var_genes) == 0 && "pca" %in% names(seurat_object@reductions)) {
          pca_loadings <- seurat_object@reductions$pca@feature.loadings
          if (!is.null(pca_loadings) && nrow(pca_loadings) > 0) {
            var_genes <- rownames(pca_loadings)
            if (verbose) message(sprintf("Extracted %d variable features from PCA loadings", length(var_genes)))
          }
        }

        if (length(var_genes) > 0) {
          var$highly_variable <- rownames(var) %in% var_genes
          # Save the order of variable features
          var$highly_variable_rank <- NA
          var$highly_variable_rank[match(var_genes, rownames(var))] <- seq_along(var_genes)
          if (verbose) message(sprintf("Marked %d variable features with order", sum(var$highly_variable)))
        }
      }, error = function(e) {
        if (verbose) message("Could not mark variable features")
      })
    }

    obs_df <- pandas$DataFrame(obs, index = rownames(obs))
    var_df <- pandas$DataFrame(var, index = rownames(var))

    # Prepare obsm (dimensional reductions) BEFORE creating AnnData
    obsm_dict <- NULL
    varm_dict <- NULL

    if (reductions && length(seurat_object@reductions) > 0) {
      if (verbose) message("Preparing dimensional reductions...")

      obsm_list <- list()
      varm_list <- list()

      for (reduction_name in names(seurat_object@reductions)) {
        tryCatch({
          reduction <- seurat_object@reductions[[reduction_name]]

          # Get embeddings
          embeddings <- tryCatch({
            Seurat::Embeddings(reduction)
          }, error = function(e) {
            reduction@cell.embeddings
          })

          if (is.null(embeddings)) {
            if (verbose) message(sprintf("  - Could not extract %s embeddings", reduction_name))
            next
          }

          # Map Seurat names to scanpy conventions
          obsm_name <- switch(reduction_name,
                              "pca" = "X_pca",
                              "umap" = "X_umap",
                              "tsne" = "X_tsne",
                              paste0("X_", reduction_name))

          # Convert to numpy array
          embeddings_py <- np$array(as.matrix(embeddings))
          obsm_list[[obsm_name]] <- embeddings_py

          if (verbose) message(sprintf("  - Prepared %s", obsm_name))

          # Save loadings for PCA (pad with zeros if needed)
          if (reduction_name == "pca") {
            tryCatch({
              loadings <- tryCatch({
                reduction@feature.loadings
              }, error = function(e) NULL)

              if (!is.null(loadings) && nrow(loadings) > 0 && ncol(loadings) > 0) {

                n_pcs <- ncol(loadings)
                n_genes_var <- nrow(var)

                if (verbose) {
                  message(sprintf("  - Processing PCA loadings: %d genes, %d PCs",
                                  nrow(loadings), n_pcs))
                }

                # Create full matrix with zeros for ALL genes
                full_loadings <- matrix(0,
                                        nrow = n_genes_var,
                                        ncol = n_pcs)
                rownames(full_loadings) <- rownames(var)
                colnames(full_loadings) <- colnames(loadings)

                # Fill in loadings for genes that were used in PCA
                matching_genes <- intersect(rownames(loadings), rownames(var))

                if (length(matching_genes) > 0) {
                  full_loadings[matching_genes, ] <- as.matrix(loadings[matching_genes, ])

                  # Convert to numpy array
                  varm_list[["PCs"]] <- np$array(full_loadings)

                  if (verbose) {
                    message(sprintf("  - Saved PCA loadings: %d/%d genes with non-zero weights",
                                    length(matching_genes), n_genes_var))
                  }
                } else {
                  if (verbose) message("  - Warning: No matching genes for PCA loadings")
                }
              }
            }, error = function(e) {
              if (verbose) message(sprintf("  - Could not save PCA loadings: %s", e$message))
            })
          }

        }, error = function(e) {
          if (verbose) warning(sprintf("Could not prepare reduction %s: %s",
                                       reduction_name, e$message))
        })
      }

      # Create dicts if we have data
      if (length(obsm_list) > 0) {
        obsm_dict <- do.call(reticulate::dict, obsm_list)
      }
      if (length(varm_list) > 0) {
        varm_dict <- do.call(reticulate::dict, varm_list)
      }
    }

    # Prepare obsp (graphs) BEFORE creating AnnData
    obsp_dict <- NULL

    if (graphs && length(seurat_object@graphs) > 0) {
      if (verbose) message("Preparing neighbor graphs...")

      obsp_list <- list()

      for (graph_name in names(seurat_object@graphs)) {
        tryCatch({
          graph <- seurat_object@graphs[[graph_name]]
          graph_t <- Matrix::t(graph)
          graph_py <- reticulate::r_to_py(as.matrix(graph_t))
          graph_sparse <- scipy_sparse$csr_matrix(graph_py)

          # Keep original graph names to preserve information
          obsp_list[[graph_name]] <- graph_sparse
          if (verbose) message(sprintf("  - Prepared %s graph", graph_name))

        }, error = function(e) {
          if (verbose) warning(sprintf("Could not prepare graph %s: %s",
                                       graph_name, e$message))
        })
      }

      if (length(obsp_list) > 0) {
        obsp_dict <- do.call(reticulate::dict, obsp_list)
      }
    }

    # Create AnnData WITH obsm, varm, obsp
    if (verbose) message("Creating AnnData object...")

    adata_args <- list(
      X = X_sparse,
      obs = obs_df,
      var = var_df
    )

    if (!is.null(obsm_dict)) {
      adata_args$obsm <- obsm_dict
    }
    if (!is.null(varm_dict)) {
      adata_args$varm <- varm_dict
    }
    if (!is.null(obsp_dict)) {
      adata_args$obsp <- obsp_dict
    }

    adata <- do.call(anndata$AnnData, adata_args)

    # Add raw counts if requested
    if (raw) {
      if (verbose) message("Adding raw counts...")

      raw_data <- tryCatch({
        SeuratObject::GetAssayData(seurat_object, assay = assay, layer = "counts")
      }, error = function(e) {
        SeuratObject::GetAssayData(seurat_object, assay = assay, slot = "counts")
      })

      if (!inherits(raw_data, "dgCMatrix")) {
        raw_data <- methods::as(raw_data, "dgCMatrix")
      }

      raw_t <- Matrix::t(raw_data)
      raw_py <- reticulate::r_to_py(as.matrix(raw_t))
      raw_sparse <- scipy_sparse$csr_matrix(raw_py)

      adata$raw <- adata
      adata$X <- raw_sparse
    }

    if (verbose) message("Writing H5AD file...")
    adata$write_h5ad(file)
    if (verbose) message("Successfully wrote ", file)

    return(file)
  }

  result <- basilisk::basiliskRun(
    env = entanglr_env,
    fun = .seurat_to_h5ad,
    seurat_object = seurat_object,
    file = normalizePath(file, mustWork = FALSE),
    assay = assay,
    X_name = X_name,
    raw = raw,
    reductions = reductions,
    graphs = graphs,
    var_features = var_features,
    scaled_data = scaled_data,
    verbose = verbose
  )

  invisible(result)
}
