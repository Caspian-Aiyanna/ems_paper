#!/usr/bin/env Rscript
# ============================================================================
# SSDM ensemble modelling for DBSCAN replicates (config-driven, resumable)
# Inputs : <project_root>/replicates_DBSCAN/<SP>_rep<k>.csv  (lon,lat)
# Env    : A-> paths.envi_after ; B-> paths.envi_before
# Outputs: results/SSDM/<RUN>/<SP>/replicates/rep<k>/  (model + rasters + metrics)
# Modes  : REPRO (cores_repro) vs FAST (cores_fast)
# NOTE   : 'rep' is a LITERAL INTEGER (10), not a variable.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(raster)    # SSDM expects Raster*
  library(SSDM)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(sf)
})

# --- robust sourcing regardless of CWD ---------------------------------------
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  filearg <- grep("^--file=", args, value = TRUE)
  if (length(filearg)) return(normalizePath(sub("^--file=", "", filearg)))
  if (!is.null(sys.frames()) && length(sys.frames())) {
    fi <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
    if (!is.na(fi)) return(fi)
  }
  stop("Cannot determine script path; run via Rscript or setwd() to project root.")
}
.script <- dirname(.this_file())
.root   <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = TRUE)

# project utils
source(file.path(.root, "R", "utils_io.R"))       # dir_ensure, log_line, %||%, safe_read_csv/safe_write_csv, skip_if_done
source(file.path(.root, "R", "utils_repro.R"))    # read_config, set_mode, seed_for
source(file.path(.root, "R", "utils_kendall.R"))  # kendall_plan

`%||%` <- function(x, y) if (is.null(x)) y else x
sp_hash_int <- function(sp) sum(utf8ToInt(as.character(sp))) %% 10000L
strict_int <- function(x, min_val = 1L) {
  xi <- suppressWarnings(as.integer(x))
  if (length(xi) != 1L || is.na(xi) || xi < min_val) min_val else xi
}
drop_nzv <- function(df, cols) {
  keep <- vapply(cols, function(cl) {
    v <- df[[cl]]
    if (is.factor(v)) return(length(unique(v)) > 1L)
    suppressWarnings(stats::sd(v, na.rm = TRUE)) > 0
  }, logical(1))
  cols[keep]
}

# --- CLI ---------------------------------------------------------------------
opt <- list(
  make_option(c("--run"),         type = "character", default = "B",   help = "Run tag: A or B"),
  make_option(c("--mode"),        type = "character", default = NULL,  help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"),     type = "character", default = NULL,  help = "Optional species filter, e.g. E5A (trains its replicates only)"),
  make_option(c("--input"),       type = "character", default = NULL,  help = "Replicates folder (default: <project_root>/replicates_DBSCAN)"),
  make_option(c("--algos"),       type = "character", default = NULL,  help = "Comma-separated algorithm list (override defaults)")
)
opts <- parse_args(OptionParser(option_list = opt))

# --- config / mode -----------------------------------------------------------
cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)

# --- paths -------------------------------------------------------------------
paths <- list(
  envi_before   = cfg$paths$envi_before   %||% file.path("data","envi","B"),
  envi_after    = cfg$paths$envi_after    %||% file.path("data","envi","A"),
  results_ssdm  = cfg$paths$results_ssdm  %||% file.path("results","SSDM"),
  plans         = cfg$paths$plans         %||% "plans",
  logs          = cfg$paths$logs          %||% "logs"
)

rep_dir <- opts$input %||% file.path(.root, "replicates_DBSCAN")

dir_ensure(paths$results_ssdm); dir_ensure(paths$plans); dir_ensure(paths$logs)
logf <- file.path(paths$logs, sprintf("04a_ssdm_replicates_%s.log", opts$run))

# --- env selection by run tag -----------------------------------------------
env_dir <- if (toupper(opts$run) == "A") paths$envi_after else paths$envi_before

log_line(sprintf("Starting 04a_ssdm_replicates.R  (mode=%s, run=%s)", mode, opts$run), logf)
log_line(sprintf("Replicates dir: %s", rep_dir), logf)
log_line(sprintf("Env dir: %s", env_dir), logf)
log_line(sprintf("Results root: %s", paths$results_ssdm), logf)

# --- load env rasters --------------------------------------------------------
ras_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(ras_files) > 0)
ras_names <- make.names(tools::file_path_sans_ext(basename(ras_files)), unique = TRUE)
Renv <- terra::rast(ras_files); names(Renv) <- ras_names

# --- Kendall keepvars --------------------------------------------------------
keepvars_csv <- file.path(paths$plans, opts$run, "keepvars.csv")
if (file.exists(keepvars_csv)) {
  keep_vars <- readr::read_csv(keepvars_csv, show_col_types = FALSE)$variable
  log_line(sprintf("Read keepvars: %d variables", length(keep_vars)), logf)
} else {
  dir_ensure(file.path(paths$plans, opts$run))
  keep_vars <- kendall_plan(opts$run, Renv,
                            out_dir = file.path(paths$plans, opts$run),
                            cutoff   = cfg$kendall$cutoff %||% 0.8,
                            nsample  = cfg$kendall$sample %||% 5000L)
  log_line(sprintf("Created keepvars: %d variables", length(keep_vars)), logf)
}
Renv <- Renv[[intersect(keep_vars, names(Renv))]]
# SSDM wants RasterStack
env_stack <- raster::stack(as(Renv, "Raster"))
names(env_stack) <- names(Renv)

# --- algorithms & cores ------------------------------------------------------
algos_default <- c("GLM","GBM","RF","MARS","SVM","ANN","CTA")
algos <- if (!is.null(opts$algos)) trimws(unlist(strsplit(opts$algos, ","))) else algos_default
algos <- intersect(algos, algos_default)
if (!length(algos)) {
  log_line("Provided --algos filtered out all defaults; falling back to full default set.", logf)
  algos <- algos_default
}

cores <- if (toupper(cfg$mode) == "REPRO") (cfg$ssdm$cores_repro %||% 1L) else (cfg$ssdm$cores_fast %||% 8L)
cores <- strict_int(cores, 1L)
log_line(sprintf("SSDM cores=%d, rep(literal)=10", cores), logf)

# --- list replicate CSVs -----------------------------------------------------
all_rep <- list.files(rep_dir, pattern = "^.+_rep[0-9]+\\.csv$", full.names = TRUE)
stopifnot(length(all_rep) > 0)

if (!is.null(opts$species)) {
  base <- paste0("^", gsub("[^A-Za-z0-9]", "", opts$species), "_rep[0-9]+\\.csv$")
  all_rep <- all_rep[grepl(base, basename(all_rep), ignore.case = FALSE)]
  if (!length(all_rep)) {
    log_line(sprintf("No replicate files matched --species %s; exiting.", opts$species), logf)
    quit(status = 0)
  }
}

# --- train per replicate -----------------------------------------------------
for (f in sort(all_rep)) {
  fname   <- basename(f)                 # e.g., "E5A_rep2.csv"
  sp      <- sub("_rep[0-9]+\\.csv$", "", fname)
  rep_idx <- as.integer(sub("^.*_rep([0-9]+)\\.csv$", "\\1", fname))

  out_dir <- file.path(paths$results_ssdm, opts$run, sp, "replicates", paste0("rep", rep_idx))
  dir_ensure(out_dir)

  # resume guard
  if (skip_if_done(out_dir, c(
    sprintf("ESDM_%s_rep%d.rds", sp, rep_idx),
    sprintf("ESDM_%s_rep%d.tif", sp, rep_idx)
    # uncertainty is optional; don't gate on it
  ))) {
    log_line(sprintf("Skip %s (rep %d) — required outputs present", sp, rep_idx), logf)
    next
  }

  log_line(sprintf("Training %s (rep %d)", sp, rep_idx), logf)

  # presences (lon, lat)
  occ_df <- readr::read_csv(f, show_col_types = FALSE)
  stopifnot(all(c("lon","lat") %in% names(occ_df)))

  # ensure predictors exist at presences
  ext_vals <- terra::extract(Renv, terra::vect(occ_df, geom = c("lon","lat"), crs = terra::crs(Renv)))
  ext_vals <- as_tibble(ext_vals)[, -1, drop = FALSE]  # drop ID column
  tmp <- dplyr::bind_cols(occ_df, as.data.frame(ext_vals)) %>% tidyr::drop_na()

  x_vars <- intersect(names(Renv), keep_vars)
  x_vars <- intersect(x_vars, names(tmp))
  x_vars <- drop_nzv(tmp, x_vars)
  if (length(x_vars) == 0L) {
    log_line(sprintf("Skip %s (rep %d) — no usable predictors at presences", sp, rep_idx), logf); next
  }

  # ensure legal, unique names
  x_vars <- make.names(x_vars, unique = TRUE)
  env_use <- raster::stack(as(Renv[[x_vars]], "Raster"))
  names(env_use) <- x_vars

  # deterministic seed per species+rep
  set.seed((cfg$seeds$per_dataset_base %||% 30000L) + sp_hash_int(paste0(sp, "_rep", rep_idx)))

  # ---- Fit ensemble (rep is LITERAL 10) -------------------------------------
  esdm <- tryCatch(
    SSDM::ensemble_modelling(
      algorithms      = algos,
      Occurrences     = tmp[, c("lon","lat")],
      Env             = env_use,
      Xcol            = "lon",
      Ycol            = "lat",
      Pcol            = NULL,              # presence-only
      rep             = 10,                # literal integer
      cores           = cores,
      cv              = "holdout",
      cv.param        = c(0.75, 10),       # literal integers
      ensemble.thresh = if (toupper(cfg$mode) == "REPRO") 0.75 else 0,
      verbose         = FALSE
    ),
    error = function(e) {
      log_line(sprintf("ensemble_modelling() error for %s (rep %d): %s",
                       sp, rep_idx, conditionMessage(e)), logf)
      NULL
    }
  )
  if (is.null(esdm)) { next }

  # save model
  saveRDS(esdm, file = file.path(out_dir, sprintf("ESDM_%s_rep%d.rds", sp, rep_idx)))

  # projection raster
  pred <- esdm@projection
  stopifnot(inherits(pred, "RasterLayer"))
  out_tif <- file.path(out_dir, sprintf("ESDM_%s_rep%d.tif", sp, rep_idx))
  raster::writeRaster(pred, filename = out_tif, format = "GTiff",
                      overwrite = TRUE, NAflag = -9999, datatype = "FLT4S")

  # uncertainty raster (if present)
  if (!is.null(esdm@uncertainty) && inherits(esdm@uncertainty, "RasterLayer")) {
    unc_tif <- file.path(out_dir, sprintf("uncertainty_%s_rep%d.tif", sp, rep_idx))
    raster::writeRaster(esdm@uncertainty, filename = unc_tif, format = "GTiff",
                        overwrite = TRUE, NAflag = -9999, datatype = "FLT4S")
  }

  # evaluation + varimp
  ev_df <- tryCatch(as.data.frame(esdm@evaluation), error = function(e) NULL)
  if (!is.null(ev_df)) {
    ev_df$dataset   <- sp
    ev_df$replicate <- rep_idx
    readr::write_csv(ev_df, file.path(out_dir, sprintf("ESDM_evaluation_%s_rep%d.csv", sp, rep_idx)))
  }

  vi <- esdm@variable.importance
  if (!is.null(vi)) {
    if (is.numeric(vi) && !is.null(names(vi))) {
      vi_df <- data.frame(Variable = names(vi), Importance = as.numeric(vi), row.names = NULL)
    } else if (is.matrix(vi) || is.data.frame(vi)) {
      vi_df <- as.data.frame(vi)
      if (!is.null(rownames(vi_df))) vi_df$Variable <- rownames(vi_df)
      if (!"Variable" %in% names(vi_df)) vi_df$Variable <- paste0("V", seq_len(nrow(vi_df)))
      if (!"Importance" %in% names(vi_df) && ncol(vi_df) >= 2) {
        names(vi_df)[2] <- "Importance"
      }
    } else vi_df <- NULL

    if (!is.null(vi_df) && nrow(vi_df) > 0) {
      readr::write_csv(vi_df, file.path(out_dir, sprintf("varimp_%s_rep%d.csv", sp, rep_idx)))
      if ("Importance" %in% names(vi_df)) {
        p <- ggplot(vi_df, aes(x = reorder(Variable, Importance), y = Importance)) +
          geom_col() + coord_flip() + theme_minimal(base_size = 11) +
          labs(title = paste("Variable importance —", sp, "(rep", rep_idx, ")"),
               x = "Variable", y = "Importance")
        ggsave(file.path(out_dir, sprintf("varimp_%s_rep%d.png", sp, rep_idx)), p,
               width = 7, height = 5, dpi = 300)
      }
    }
  }

  log_line(sprintf("Done %s (rep %d)", sp, rep_idx), logf)
}

log_line("All replicate SSDM jobs finished.", logf)
