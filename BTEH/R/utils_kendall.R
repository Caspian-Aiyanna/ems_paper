suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)   # we call tibble() below
})

# Greedy Kendall-pruning: drop vars until all |tau| < cutoff
prune_by_kendall <- function(cmat, cutoff = 0.8) {
  M <- abs(cmat)
  diag(M) <- 0
  if (ncol(M) < 2) return(character(0))
  to_drop <- character(0)
  keep    <- colnames(M)
  repeat {
    maxv <- suppressWarnings(max(M, na.rm = TRUE))
    if (!is.finite(maxv) || maxv < cutoff) break
    idx  <- which(M == maxv, arr.ind = TRUE)[1, ]
    c1   <- colnames(M)[idx[1]]
    c2   <- colnames(M)[idx[2]]
    m1 <- suppressWarnings(mean(M[idx[1], -idx[1]], na.rm = TRUE))
    m2 <- suppressWarnings(mean(M[idx[2], -idx[2]], na.rm = TRUE))
    drop_var <- if (is.na(m1) || is.na(m2)) c2 else if (m1 >= m2) c1 else c2
    to_drop  <- c(to_drop, drop_var)
    keep <- setdiff(keep, drop_var)
    M <- M[keep, keep, drop = FALSE]
    if (length(keep) < 2) break
  }
  unique(to_drop)
}

# Deterministic cell sampling helper (not strictly required but handy)
sample_cells_deterministic <- function(r, n, seed = 123L) {
  set.seed(seed)
  r1 <- r[[1]]
  pts <- spatSample(r1, size = min(n, ncell(r1)), method = "random", as.points = TRUE, na.rm = TRUE)
  crds(pts)
}

# Main: compute Kendall, write plans/<run>/keepvars.csv etc., return kept vars
kendall_plan <- function(run,
                         r,
                         out_dir = file.path("plans", run),
                         cutoff = 0.8,
                         nsample = 5000L) {
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  nsamp <- min(as.integer(nsample), ncell(r[[1]]))
  pts   <- spatSample(r[[1]], size = nsamp, method = "random", as.points = TRUE, na.rm = TRUE)
  vals  <- terra::extract(r, pts)[, -1, drop = FALSE]

  df <- as.data.frame(vals)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) > nsample) {
    set.seed(1L)
    df <- df[sample(nrow(df), nsample), , drop = FALSE]
  }

  kt <- suppressWarnings(stats::cor(df, method = "kendall", use = "pairwise.complete.obs"))
  drop_vars <- prune_by_kendall(kt, cutoff = cutoff)
  keep_vars <- setdiff(colnames(kt), drop_vars)

  readr::write_csv(tibble(variable = keep_vars), file.path(out_dir, "keepvars.csv"))
  readr::write_csv(tibble(dropped_variable = drop_vars), file.path(out_dir, "kendall_dropped_vars.csv"))

  # Heatmap of kept vars
  if (length(keep_vars) >= 1) {
    kt_keep <- kt[keep_vars, keep_vars, drop = FALSE]
    kt_long <- as.data.frame(as.table(kt_keep))
    names(kt_long) <- c("Var1", "Var2", "tau")
    ggsave(file.path(out_dir, "kendall_heatmap_kept.png"),
           ggplot(kt_long, aes(Var1, Var2, fill = tau)) +
             geom_tile() +
             scale_fill_gradient2(limits = c(-1, 1)) +
             labs(x = NULL, y = NULL, fill = "Kendall τ",
                  title = paste0("Kendall (|τ|<", cutoff, ") — ", run)) +
             theme_minimal(base_size = 10) +
             theme(axis.text.x = element_text(angle = 90, hjust = 1)),
           width = 8, height = 6, dpi = 300)
  }

  return(keep_vars)
}
