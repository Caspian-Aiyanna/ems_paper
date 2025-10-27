# ============================================================================
# Purpose: H2O helpers (deterministic vs fast modes, nthreads incl. "-2" logic),
#          prediction over rasters in chunks, and simple IO wrappers.
# ============================================================================

suppressPackageStartupMessages({
  library(h2o)
  library(terra)
})

# Reserve cores logic: if nthreads is negative, reserve abs(nthreads) cores.
# e.g., -2 -> use (available_cores - 2), minimum 1.
resolve_nthreads <- function(nthreads_cfg) {
  if (is.null(nthreads_cfg)) return(1L)
  if (nthreads_cfg > 0) return(as.integer(nthreads_cfg))
  # negative -> reserve cores
  avail <- as.integer(parallel::detectCores(logical = TRUE))
  use   <- max(1L, avail - abs(as.integer(nthreads_cfg)))
  return(use)
}

init_h2o <- function(cfg) {
  mode <- toupper(cfg$mode %||% "REPRO")
  if (mode == "REPRO") {
    nthreads <- resolve_nthreads(cfg$h2o$nthreads_repro %||% 1L)
  } else {
    nthreads <- resolve_nthreads(cfg$h2o$nthreads_fast %||% -2L)
  }

  # Quietly shut down any existing cluster (avoids version mismatches)
  suppressWarnings(try(h2o::h2o.shutdown(prompt = FALSE), silent = TRUE))
  Sys.sleep(1)

  # Start fresh cluster; no global seed API in 3.46.x (seed is passed to algos)
  h2o::h2o.init(nthreads = nthreads, strict_version_check = FALSE)
  h2o::h2o.removeAll()

  invisible(nthreads)
}



# Predict a raster stack using an H2O binomial model.
# Processes in blocks to manage memory. Returns a SpatRaster.
# Assumes raster names exactly match training feature names order.
predict_raster_h2o <- function(r, model, block_rows = 50000) {
  stopifnot(inherits(r, "SpatRaster"))
  # prepare an empty output raster with same geometry
  out <- rast(r[[1]])
  names(out) <- "p1"

  ncell_total <- ncell(r)
  idx <- seq_len(ncell_total)

  # terra::as.data.frame with cells argument to avoid huge memory spikes
  # We'll iterate by chunks of row indices in the cell vector.
  vals <- rep(NA_real_, ncell_total)

  # Chunking helper
  chunks <- split(idx, ceiling(seq_along(idx) / block_rows))
  for (ch in chunks) {
    df <- as.data.frame(r[ch], na.rm = FALSE)
    hf <- as.h2o(df)
    pr <- h2o.predict(model, hf)
    # 'p1' column for binomial
    vals[ch] <- as.vector(pr[["p1"]])
  }
  out <- setValues(out, vals)
  return(out)
}

# Simple guard to skip heavy work if all files already exist
skip_if_done <- function(dir, files) {
  all(file.exists(file.path(dir, files)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a