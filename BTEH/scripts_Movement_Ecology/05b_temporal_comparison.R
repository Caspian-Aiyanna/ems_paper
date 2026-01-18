# =============================================================================
# 05b_temporal_comparison.R — Temporal Change Analysis (Before vs After)
#
# Description: This script quantifies changes in elephant habitat suitability
# before and after fence removal. It calculates area metrics, Jaccard overlap,
# and area gain/loss across various suitability thresholds.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(sf)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
dir_ensure <- function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
log_line <- function(msg, file = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, msg))
  if (!is.null(file)) write(paste0("[", ts, "] ", msg), file, append = TRUE)
}

# ---------- robust script path (for self-invocation) --------------------------
.this_file <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  if (length(f)) {
    return(normalizePath(f, mustWork = FALSE))
  }
  if (!is.null(sys.frames()) && length(sys.frames())) {
    fi <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (!is.null(fi)) {
      return(normalizePath(fi, mustWork = FALSE))
    }
  }
  NA_character_
}
.script <- dirname(.this_file())
if (!nzchar(.script)) .script <- "."
.root <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = FALSE)

# ---------- CLI ---------------------------------------------------------------
opt_list <- list(
  make_option(c("--run"),
    type = "character", default = NULL,
    help = "Run tag: A or B. If omitted, the script runs BOTH A and B."
  ),
  make_option(c("--ids"),
    type = "character", default = "E3,E4,E5",
    help = "Comma-separated base IDs (no A/B), e.g. 'E3,E4,E5'."
  ),
  make_option(c("--t"),
    type = "character", default = "0.25,0.50,0.75",
    help = "Comma-separated absolute thresholds on 0–1 scale."
  ),
  make_option(c("--aoi"),
    type = "character", default = "data/shp/HV20233.shp",
    help = "Optional AOI vector for plots; set '' to disable."
  ),
  make_option(c("--plots"),
    type = "character", default = "TRUE",
    help = "TRUE/FALSE: write summary plots."
  )
)
opts <- parse_args(OptionParser(option_list = opt_list))

# ---------- If --run omitted, run BOTH sequentially --------------------------
if (is.null(opts$run)) {
  script <- .this_file()
  if (!nzchar(script)) stop("[FATAL] Cannot resolve script path for self-invocation.")
  args_in <- commandArgs(trailingOnly = TRUE)
  # drop any --run/-r if present (defensive)
  drop <- which(args_in %in% c("--run", "-r"))
  if (length(drop)) args_in <- args_in[setdiff(seq_along(args_in), c(drop, drop + 1L))]
  for (r in c("A", "B")) {
    log_line(sprintf("[INFO] Launching child for run %s ...", r))
    system2(file.path(R.home("bin"), "Rscript"),
      args = c(script, args_in, "--run", r)
    )
  }
  quit(status = 0)
}

# ---------- Parse/normalize options ------------------------------------------
run_tag <- toupper(opts$run) # A or B (used only for output folder)
base_ids <- toupper(trimws(unlist(strsplit(opts$ids, "[,;\\s]+"))))
thresholds <- as.numeric(trimws(unlist(strsplit(opts$t, "[,;\\s]+"))))
make_plots <- toupper(opts$plots) %in% c("TRUE", "T", "1", "YES", "Y")

# ---------- Fixed roots -------------------------------------------------------
h2o_root <- file.path(.root, "results", "H2O")
ssdm_root <- file.path(.root, "results", "SSDM")
compare_root <- file.path(.root, "results", "compare", "h2o_vs_ssdm", run_tag)

temporal_out <- file.path(compare_root, "10_temporal")
dir_ensure(temporal_out)
logf <- file.path(temporal_out, sprintf("temporal_change_%s.log", run_tag))
log_line(sprintf("Temporal change script start (run=%s)", run_tag), logf)

# ---------- Optional AOI ------------------------------------------------------
aoi <- NULL
if (!is.null(opts$aoi) && nzchar(opts$aoi) && file.exists(opts$aoi)) {
  aoi <- tryCatch(suppressWarnings(sf::st_read(opts$aoi, quiet = TRUE)), error = function(e) NULL)
  if (!is.null(aoi)) log_line(sprintf("[AOI] Loaded %s", opts$aoi), logf)
}

# ---------- Geospatial helpers ------------------------------------------------
align_to <- function(r1, r2, categorical = FALSE) {
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- project(r2, crs(r1), method = if (categorical) "near" else "bilinear")
    r2 <- resample(r2, r1, method = if (categorical) "near" else "bilinear")
    r2 <- crop(r2, r1)
  }
  m <- !is.na(r1) & !is.na(r2)
  r1 <- mask(r1, m, maskvalues = 0)
  r2 <- mask(r2, m, maskvalues = 0)
  list(r1 = r1, r2 = r2)
}
binary_area_km2 <- function(bin_r) {
  cs <- terra::cellSize(bin_r, unit = "km")
  as.numeric(global(cs * (bin_r == 1), "sum", na.rm = TRUE)[[1]] %||% 0)
}

# ---------- I/O helpers -------------------------------------------------------
load_pair <- function(method_root, base_id, file_stub, runB = "B", runA = "A") {
  spB <- paste0(base_id, runB)
  spA <- paste0(base_id, runA)
  fB <- file.path(method_root, runB, spB, sprintf(file_stub, spB))
  fA <- file.path(method_root, runA, spA, sprintf(file_stub, spA))
  if (!file.exists(fB) || !file.exists(fA)) {
    return(NULL)
  }
  list(B = rast(fB), A = rast(fA), spB = spB, spA = spA)
}

temporal_metrics_once <- function(rB, rA, cutoff) {
  al <- align_to(rB, rA, categorical = FALSE)
  rB <- al$r1
  rA <- al$r2
  bB <- rB > cutoff
  bA <- rA > cutoff

  loss <- (bB == 1) & (bA == 0)
  stable <- (bB == 1) & (bA == 1)
  gain <- (bB == 0) & (bA == 1)

  B_area <- binary_area_km2(bB)
  A_area <- binary_area_km2(bA)
  Loss <- binary_area_km2(loss)
  Stable <- binary_area_km2(stable)
  Gain <- binary_area_km2(gain)

  pct_delta <- if (B_area > 0) (A_area - B_area) / B_area * 100 else NA_real_
  inter <- as.numeric(global((bB == 1) & (bA == 1), "sum", na.rm = TRUE)[[1]])
  union <- as.numeric(global((bB == 1) | (bA == 1), "sum", na.rm = TRUE)[[1]])
  Jt <- ifelse(union == 0, NA_real_, inter / union)

  tibble(
    B_area_km2 = B_area, A_area_km2 = A_area,
    Loss_km2 = Loss, Stable_km2 = Stable, Gain_km2 = Gain,
    pct_delta = pct_delta, J_t = Jt
  )
}

# ---------- Main: loop IDs and methods ---------------------------------------
rows <- list()
for (base in base_ids) {
  # SSDM
  ssdm_pair <- load_pair(ssdm_root, base, "ESDM_%s.tif")
  if (is.null(ssdm_pair)) log_line(sprintf("[WARN] Missing SSDM pair for %s", base), logf)
  if (!is.null(ssdm_pair)) {
    for (tcut in thresholds) {
      met <- temporal_metrics_once(ssdm_pair$B, ssdm_pair$A, tcut)
      rows[[length(rows) + 1L]] <- mutate(met, ID = base, Method = "SSDM", t = tcut, .before = 1)
    }
  }

  # H2O
  h2o_pair <- load_pair(h2o_root, base, "prediction_%s.tif")
  if (is.null(h2o_pair)) log_line(sprintf("[WARN] Missing H2O pair for %s", base), logf)
  if (!is.null(h2o_pair)) {
    for (tcut in thresholds) {
      met <- temporal_metrics_once(h2o_pair$B, h2o_pair$A, tcut)
      rows[[length(rows) + 1L]] <- mutate(met, ID = base, Method = "h2o", t = tcut, .before = 1)
    }
  }
}

temporal_df <- bind_rows(rows) %>%
  select(
    ID, Method, t,
    B_area_km2, A_area_km2,
    Loss_km2, Stable_km2, Gain_km2,
    pct_delta, J_t
  )

# ---------- Write CSV ---------------------------------------------------------
csv_out <- file.path(temporal_out, "temporal_change_metrics.csv")
readr::write_csv(temporal_df, csv_out)
log_line(sprintf("Wrote CSV: %s", csv_out), logf)

# ---------- Quick plots (optional) -------------------------------------------
if (make_plots && nrow(temporal_df)) {
  p_j <- ggplot(temporal_df, aes(x = factor(t), y = J_t, group = Method, linetype = Method, color = Method)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ID, nrow = 1) +
    labs(x = "Threshold (t)", y = "Jaccard(A,B)", title = "Temporal overlap across thresholds") +
    theme_minimal(base_size = 11)
  ggsave(file.path(temporal_out, "jaccard_vs_threshold.png"), p_j, width = 9.5, height = 3.4, dpi = 300)

  p_pc <- ggplot(temporal_df, aes(x = factor(t), y = pct_delta, fill = Method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    facet_wrap(~ID, nrow = 1) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(x = "Threshold (t)", y = "%Δ Area (A vs B)", title = "Temporal area change by threshold") +
    theme_minimal(base_size = 11)
  ggsave(file.path(temporal_out, "percent_area_change_by_threshold.png"), p_pc, width = 9.5, height = 3.4, dpi = 300)
}

log_line(sprintf("Done. Outputs → %s", temporal_out), logf)
