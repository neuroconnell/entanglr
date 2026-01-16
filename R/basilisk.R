#' @importFrom basilisk BasiliskEnvironment basiliskRun
NULL

# Python environment for entanglr
entanglr_env <- basilisk::BasiliskEnvironment(
  envname = "entanglr_env",
  pkgname = "entanglr",
  packages = c("python==3.9"),
  pip = c(
    "anndata==0.9.2",
    "scanpy==1.9.3",
    "numpy==1.23.5",
    "pandas==1.5.3",
    "h5py==3.8.0"
  )
)
