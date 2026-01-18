# =============================================================================
# 04_ssdm_train.R — Stacked Species Distribution Modeling (SSDM) Pipeline
#
# Description: This script trains ensemble species distribution models using the
# SSDM package. It integrates multiple algorithms (GLM, GBM, RF, etc.) and
# generates suitability maps and evaluation metrics for each elephant.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(raster) # SSDM expects Raster*
  library(SSDM)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(pROC)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
sp_hash_int <- function(sp) sum(utf8ToInt(as.character(sp))) %% 10000L
drop_nzv <- function(df, cols) {
  keep <- vapply(cols, function(cl) {
    v <- df[[cl]]
    if (is.factor(v)) {
      return(length(unique(v)) > 1L)
    }
    suppressWarnings(stats::sd(v, na.rm = TRUE)) > 0
  }, logical(1))
  cols[keep]
}
strict_int <- function(x, min_val = 1L) {
  xi <- suppressWarnings(as.integer(x))
  if (length(xi) != 1L || is.na(xi) || xi < min_val) min_val else xi
}
nz <- function(x) {
  x <- if (is.null(x)) NA_character_ else as.character(x)
  x[!nzchar(x)] <- NA_character_
  x
}

# --- robust sourcing regardless of CWD ---------------------------------------
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  filearg <- grep("^--file=", args, value = TRUE)
  if (length(filearg)) {
    return(normalizePath(sub("^--file=", "", filearg)))
  }
  if (!is.null(sys.frames()) && length(sys.frames())) {
    fi <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
    if (!is.na(fi)) {
      return(fi)
    }
  }
  stop("Cannot determine script path; run via Rscript or setwd() to project root.")
}
.script <- dirname(.this_file())
.root <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = TRUE)

source(file.path(.root, "R", "utils_io.R")) # dir_ensure, log_line, %||%
source(file.path(.root, "R", "utils_repro.R")) # read_config, set_mode
source(file.path(.root, "R", "utils_kendall.R")) # kendall_plan (for keepvars)

# --- CLI ---------------------------------------------------------------------
opt <- list(
  # Set default to "A" or "B" to bake-in your preferred default;
  # leave NULL to force providing --run or using ENV / cfg fallback.
  make_option(c("--run"), type = "character", default = "A", help = "Run tag: A or B."),
  make_option(c("--mode"), type = "character", default = "REPRO", help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"), type = "character", default = NULL, help = "Optional single species (e.g., E3A)"),
  make_option(c("--algos"), type = "character", default = NULL, help = "Comma-separated algorithms (override)"),
  make_option(c("--spatial_cv"), action = "store_true", default = FALSE, help = "Ignored (holdout only)")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg) # "REPRO" or "FAST"

# --- resolve run (CLI > Option default > ENV > cfg; else error) ---------------
resolve_run <- function(opt_val, cfg) {
  cand <- c(
    nz(opt_val)[1],
    nz(Sys.getenv("RUN", unset = NA_character_))[1],
    nz(cfg$defaults$run %||% NA_character_)[1]
  )
  cand <- toupper(trimws(cand))
  cand <- cand[!is.na(cand)]
  if (length(cand) == 0L) {
    stop(
      "No run specified. Choose A or B via one of:\n",
      " - pass --run A|B\n",
      " - set Option default to 'A' or 'B' in this script\n",
      " - set environment variable RUN=A|B\n",
      " - set cfg$defaults$run: A|B"
    )
  }
  run <- cand[1]
  if (!run %in% c("A", "B")) stop("`run` must be A or B (got: ", run, ")")
  run
}
run_tag <- resolve_run(opts$run, cfg)
message(sprintf("[run=%s] Using run determined via precedence (CLI > Option default > ENV > cfg).", run_tag))

# --- paths & logs ------------------------------------------------------------
occ_dir <- cfg$paths$occ %||% file.path(.root, "data", "occ", "thinned_DBSCAN")
res_root <- cfg$paths$results_ssdm %||% file.path(.root, "results", "SSDM")
plans_dir <- cfg$paths$plans %||% file.path(.root, "plans")
logs_dir <- cfg$paths$logs %||% file.path(.root, "logs")
dir_ensure(res_root)
dir_ensure(plans_dir)
dir_ensure(logs_dir)

logf <- file.path(logs_dir, sprintf("04_ssdm_train_%s.log", run_tag))
log_line(sprintf("Starting 04_ssdm_train.R  (mode=%s, run=%s)", mode, run_tag), logf)

# --- choose raster dir by run tag (explicit, no fallback) --------------------
env_dir <- if (run_tag == "A") {
  cfg$paths$envi_after %||% file.path("data", "envi", "A")
} else {
  cfg$paths$envi_before %||% file.path("data", "envi", "B")
}
if (!dir.exists(env_dir)) stop("Env dir not found for run=", run_tag, ": ", env_dir)
ras_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(ras_files) == 0L) stop("No .tif in ", env_dir, " (run=", run_tag, ")")
ras_names <- make.names(tools::file_path_sans_ext(basename(ras_files)), unique = TRUE)
Renv <- terra::rast(ras_files)
names(Renv) <- ras_names

# --- Kendall keepvars --------------------------------------------------------
keepvars_csv <- file.path(plans_dir, run_tag, "keepvars.csv")
if (file.exists(keepvars_csv)) {
  keep_vars <- readr::read_csv(keepvars_csv, show_col_types = FALSE)$variable
  log_line(sprintf("Read keepvars: %d variables", length(keep_vars)), logf)
} else {
  dir_ensure(file.path(plans_dir, run_tag))
  keep_vars <- kendall_plan(
    run_tag, Renv,
    out_dir = file.path(plans_dir, run_tag),
    cutoff = cfg$kendall$cutoff %||% 0.8,
    nsample = cfg$kendall$sample %||% 5000L
  )
  log_line(sprintf("Created keepvars: %d variables", length(keep_vars)), logf)
}
K <- intersect(keep_vars, names(Renv))
stopifnot(length(K) > 0)
Renv <- Renv[[K]]

# --- cores & replications (resolve once) -------------------------------------
cores_repro <- strict_int(cfg$ssdm$cores_repro %||% 1L, 1L)
cores_fast <- strict_int(cfg$ssdm$cores_fast %||% 8L, 1L)
replications <- 10L # literal per your constraint

log_line(sprintf("SSDM cores: REPRO=%d, FAST=%d; replications=%d", cores_repro, cores_fast, replications), logf)

# --- species list ------------------------------------------------------------
occ_files <- list.files(occ_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(occ_files) > 0)
sp_list <- sort(tools::file_path_sans_ext(basename(occ_files)))
if (!is.null(opts$species)) {
  sp_list <- intersect(sp_list, opts$species)
  if (length(sp_list) == 0L) stop("Requested --species not found in ", occ_dir, ": ", opts$species)
}
log_line(sprintf("Species to run (%d): %s", length(sp_list), paste(sp_list, collapse = ", ")), logf)

# --- algorithms --------------------------------------------------------------
algos_default <- c("GLM", "GBM", "RF", "MARS", "SVM", "ANN", "CTA")
algos <- if (!is.null(opts$algos)) trimws(unlist(strsplit(opts$algos, ","))) else algos_default

# --- ESDM trainers (use resolved cores; literal rep/cv/thresholds) -----------
train_esdm_repro <- function(tag, occ_df, env_use, out_dir, algos_eff, cores_eff) {
  dir_ensure(out_dir)
  message(sprintf(
    "ESDM training [REPRO]: %s | n=%d | vars=%d | cores=%d | reps=%d",
    tag, nrow(occ_df), raster::nlayers(env_use), cores_eff, 10
  ))
  log_line(sprintf("Training ESDM %s (REPRO: cores=%d, rep=10)", tag, cores_eff), logf)

  occ_use <- occ_df %>%
    dplyr::select(lon, lat) %>%
    dplyr::distinct()
  set.seed((cfg$seeds$per_dataset_base %||% 30000L) + sp_hash_int(tag))

  esdm <- tryCatch(
    SSDM::ensemble_modelling(
      algorithms      = algos_eff,
      Occurrences     = occ_use,
      Env             = env_use,
      Xcol            = "lon",
      Ycol            = "lat",
      Pcol            = NULL,
      rep             = 10,
      cores           = cores_eff,
      cv              = "holdout",
      cv.param        = c(0.75, 10),
      ensemble.thresh = 0.75,
      verbose         = FALSE
    ),
    error = function(e) {
      log_line(sprintf("ensemble_modelling() error for %s: %s", tag, conditionMessage(e)), logf)
      NULL
    }
  )
  if (is.null(esdm)) stop(sprintf("ESDM training failed for %s", tag))
  save(esdm, file = file.path(out_dir, sprintf("ESDM_%s.rda", tag)))
  raster::writeRaster(esdm@projection,
    filename = file.path(out_dir, sprintf("ESDM_%s.tif", tag)),
    format = "GTiff", overwrite = TRUE
  )

  ev_df <- tryCatch(as.data.frame(esdm@evaluation), error = function(e) NULL)
  if (!is.null(ev_df)) {
    ev_df$dataset <- tag
    readr::write_csv(ev_df, file.path(out_dir, sprintf("ESDM_evaluation_%s.csv", tag)))
  }

  vi <- esdm@variable.importance
  if (!is.null(vi)) {
    if (is.numeric(vi) && !is.null(names(vi))) {
      vi_df <- data.frame(Variable = names(vi), Importance = as.numeric(vi), row.names = NULL)
    } else if (is.matrix(vi) || is.data.frame(vi)) {
      tmp <- as.data.frame(vi)
      if (!is.null(rownames(tmp))) tmp$Variable <- rownames(tmp)
      if (!"Variable" %in% names(tmp)) tmp$Variable <- paste0("V", seq_len(nrow(tmp)))
      long_vi <- tidyr::pivot_longer(tmp, cols = setdiff(names(tmp), "Variable"), names_to = "Model", values_to = "Importance")
      vi_df <- long_vi %>%
        dplyr::group_by(Variable) %>%
        dplyr::summarise(Importance = mean(Importance, na.rm = TRUE), .groups = "drop")
    } else {
      vi_df <- NULL
    }
    if (!is.null(vi_df) && nrow(vi_df) > 0) {
      readr::write_csv(vi_df, file.path(out_dir, sprintf("varimp_%s.csv", tag)))
      p <- ggplot(vi_df, aes(x = reorder(Variable, Importance), y = Importance)) +
        geom_col() +
        coord_flip() +
        theme_minimal(base_size = 11) +
        labs(title = paste("Variable importance —", tag), x = "Variable", y = "Importance")
      ggsave(file.path(out_dir, sprintf("varimp_%s.png", tag)), p, width = 7, height = 5, dpi = 300)
    }
  }
  invisible(esdm)
}

train_esdm_fast <- function(tag, occ_df, env_use, out_dir, algos_eff, cores_eff) {
  dir_ensure(out_dir)
  message(sprintf(
    "ESDM training [FAST]: %s | n=%d | vars=%d | cores=%d | reps=%d",
    tag, nrow(occ_df), raster::nlayers(env_use), cores_eff, 1
  ))
  log_line(sprintf("Training ESDM %s (FAST: cores=%d, rep=10)", tag, cores_eff), logf)

  occ_use <- occ_df %>%
    dplyr::select(lon, lat) %>%
    dplyr::distinct()
  set.seed((cfg$seeds$per_dataset_base %||% 30000L) + sp_hash_int(tag))

  esdm <- tryCatch(
    SSDM::ensemble_modelling(
      algorithms      = algos_eff,
      Occurrences     = occ_use,
      Env             = env_use,
      Xcol            = "lon",
      Ycol            = "lat",
      Pcol            = NULL,
      rep             = 1,
      cores           = cores_eff,
      cv              = "holdout",
      cv.param        = c(0.75, 10),
      ensemble.thresh = 0,
      verbose         = FALSE
    ),
    error = function(e) {
      log_line(sprintf("ensemble_modelling() error for %s: %s", tag, conditionMessage(e)), logf)
      NULL
    }
  )
  if (is.null(esdm)) stop(sprintf("ESDM training failed for %s", tag))
  save(esdm, file = file.path(out_dir, sprintf("ESDM_%s.rda", tag)))
  raster::writeRaster(esdm@projection,
    filename = file.path(out_dir, sprintf("ESDM_%s.tif", tag)),
    format = "GTiff", overwrite = TRUE
  )

  ev_df <- tryCatch(as.data.frame(esdm@evaluation), error = function(e) NULL)
  if (!is.null(ev_df)) {
    ev_df$dataset <- tag
    readr::write_csv(ev_df, file.path(out_dir, sprintf("ESDM_evaluation_%s.csv", tag)))
  }

  vi <- esdm@variable.importance
  if (!is.null(vi)) {
    if (is.numeric(vi) && !is.null(names(vi))) {
      vi_df <- data.frame(Variable = names(vi), Importance = as.numeric(vi), row.names = NULL)
    } else if (is.matrix(vi) || is.data.frame(vi)) {
      tmp <- as.data.frame(vi)
      if (!is.null(rownames(tmp))) tmp$Variable <- rownames(tmp)
      if (!"Variable" %in% names(tmp)) tmp$Variable <- paste0("V", seq_len(nrow(tmp)))
      long_vi <- tidyr::pivot_longer(tmp, cols = setdiff(names(tmp), "Variable"), names_to = "Model", values_to = "Importance")
      vi_df <- long_vi %>%
        dplyr::group_by(Variable) %>%
        dplyr::summarise(Importance = mean(Importance, na.rm = TRUE), .groups = "drop")
    } else {
      vi_df <- NULL
    }
    if (!is.null(vi_df) && nrow(vi_df) > 0) {
      readr::write_csv(vi_df, file.path(out_dir, sprintf("varimp_%s.csv", tag)))
      p <- ggplot(vi_df, aes(x = reorder(Variable, Importance), y = Importance)) +
        geom_col() +
        coord_flip() +
        theme_minimal(base_size = 11) +
        labs(title = paste("Variable importance —", tag), x = "Variable", y = "Importance")
      ggsave(file.path(out_dir, sprintf("varimp_%s.png", tag)), p, width = 7, height = 5, dpi = 300)
    }
  }
  invisible(esdm)
}

# --- main loop ---------------------------------------------------------------
for (sp in sp_list) {
  sp_csv <- file.path(occ_dir, paste0(sp, ".csv"))
  if (!file.exists(sp_csv)) next

  out_dir <- file.path(res_root, run_tag, sp) # keep A/B separated like H2O
  dir_ensure(out_dir)

  # resume guard
  if (skip_if_done(out_dir, c(sprintf("ESDM_%s.rda", sp), sprintf("ESDM_%s.tif", sp)))) {
    log_line(sprintf("Skip %s — outputs present", sp), logf)
    next
  }

  log_line(sprintf("Training %s", sp), logf)

  occ_df <- readr::read_csv(sp_csv, show_col_types = FALSE)
  stopifnot(all(c("lon", "lat") %in% names(occ_df)))

  # extract env values at presences to drop NA rows
  vals <- terra::extract(Renv, terra::vect(occ_df, geom = c("lon", "lat"), crs = terra::crs(Renv)))[, -1, drop = FALSE]
  tmp <- dplyr::bind_cols(occ_df, as.data.frame(vals)) %>% tidyr::drop_na()

  # predictors = keepvars present at presences; then drop zero-variance-at-presences
  x_vars <- intersect(names(Renv), keep_vars)
  x_vars <- intersect(x_vars, names(tmp))
  x_vars <- drop_nzv(tmp, x_vars)
  if (length(x_vars) == 0L) {
    log_line(sprintf("Skip %s — no usable predictors at presence points", sp), logf)
    next
  }

  # build RasterStack from aligned terra stack
  Renv_use <- Renv[[x_vars]]
  env_use <- raster::stack(as(Renv_use, "Raster"))
  names(env_use) <- names(Renv_use)

  # Optional, dataset-specific tweaks (kept, but rep stays 10 by spec)
  if (identical(sp, "E5A")) {
    cores_eff <- min(2L, if (mode == "REPRO") cores_repro else cores_fast)
    algos_eff <- intersect(algos, c("GLM", "GBM", "RF"))
  } else {
    cores_eff <- if (mode == "REPRO") cores_repro else cores_fast
    algos_eff <- algos
  }

  # Train
  if (mode == "REPRO") {
    esdm <- train_esdm_repro(
      tag = sp, occ_df = tmp[, c("lon", "lat")],
      env_use = env_use, out_dir = out_dir,
      algos_eff = algos_eff, cores_eff = cores_eff
    )
  } else {
    esdm <- train_esdm_fast(
      tag = sp, occ_df = tmp[, c("lon", "lat")],
      env_use = env_use, out_dir = out_dir,
      algos_eff = algos_eff, cores_eff = cores_eff
    )
  }

  log_line(sprintf("Done %s", sp), logf)
}

log_line("All species finished.", logf)
