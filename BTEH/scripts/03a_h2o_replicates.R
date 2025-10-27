#!/usr/bin/env Rscript
# ============================================================================
# Purpose: Train H2O models for DBSCAN replicate occurrence files
#          Mirrors 03_h2o_train.R logic (config-driven, REPRO/FAST modes,
#          Kendall keepvars, deterministic seeds, resume guards).
#          Inputs: <project_root>/replicates_DBSCAN/<SP>_rep<k>.csv
#          Outputs: results/H2O/<RUN>/<SP>/replicates/rep<k>/
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(h2o)
  library(dplyr)
  library(readr)
  library(sf)
  library(tidyr)
})

# --- robust sourcing regardless of CWD ---
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

source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))
source(file.path(.root, "R", "utils_kendall.R"))
source(file.path(.root, "R", "utils_h2o.R"))

# --- CLI ---
opt <- list(
  make_option(c("--run"),              type = "character", default = "B",   help = "Run tag: A or B (chooses envi_after/envi_before)"),
  make_option(c("--mode"),             type = "character", default = NULL,  help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"),          type = "character", default = NULL,  help = "Optional species prefix (e.g., E5A). If set, only its replicates are trained"),
  make_option(c("--input"),            type = "character", default = NULL,  help = "Folder of replicate CSVs (default: <project_root>/replicates_DBSCAN)"),
  make_option(c("--max_models"),       type = "integer",   default = NULL,  help = "Override AutoML max_models"),
  make_option(c("--max_runtime_secs"), type = "integer",   default = NULL,  help = "Optional AutoML time budget (seconds)")
)
opts <- parse_args(OptionParser(option_list = opt))

# --- config / mode ---
cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)

# --- paths (follow your config keys exactly) ---
paths <- list(
  root          = cfg$paths$root          %||% ".",
  envi_before   = cfg$paths$envi_before   %||% file.path("data","envi","B"),
  envi_after    = cfg$paths$envi_after    %||% file.path("data","envi","A"),
  occ_main      = cfg$paths$occ           %||% file.path("data","occ","thinned_DBSCAN"),
  results_h2o   = cfg$paths$results_h2o   %||% file.path("results","H2O"),
  results_ssdm  = cfg$paths$results_ssdm  %||% file.path("results","SSDM"),
  plans         = cfg$paths$plans         %||% "plans",
  clean         = cfg$paths$clean         %||% file.path("data","clean"),
  raw           = cfg$paths$raw           %||% file.path("data","raw"),
  logs          = cfg$paths$logs          %||% "logs"
)

# Replicates directory default (created by 02_dbscan_thin_degrees.R):
rep_dir <- opts$input %||% file.path(.root, "replicates_DBSCAN")

logs_dir <- file.path(paths$logs)
dir_ensure(logs_dir)
logf <- file.path(logs_dir, sprintf("03a_h2o_replicates_%s.log", opts$run))

# --- env dir chosen by RUN tag (same convention as 03_h2o_train.R) ---
env_dir <- if (toupper(opts$run) == "A") paths$envi_after else paths$envi_before

# --- banner ---
log_line(sprintf("Starting 03a_h2o_replicates.R  (mode=%s, run=%s)", mode, opts$run), logf)
log_line(sprintf("Replicates dir: %s", rep_dir), logf)
log_line(sprintf("Env dir: %s", env_dir), logf)
log_line(sprintf("Results root: %s", paths$results_h2o), logf)

# --- load environment rasters ---
ras_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(ras_files) > 0)
ras_names <- make.names(tools::file_path_sans_ext(basename(ras_files)), unique = TRUE)
Renv <- terra::rast(ras_files); names(Renv) <- ras_names

# --- Kendall plan (per RUN, shared with 03) ---
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
Renv <- Renv[[keep_vars]]

# --- H2O init ---
threads <- init_h2o(cfg)
log_line(sprintf("H2O threads: %d", threads), logf)

# --- collect replicate CSVs ---
all_rep <- list.files(rep_dir, pattern = "^.+_rep[0-9]+\\.csv$", full.names = TRUE)
stopifnot(length(all_rep) > 0)

# Filter by species if requested (match prefix before _repN)
if (!is.null(opts$species)) {
  base <- paste0("^", gsub("[^A-Za-z0-9]", "", opts$species), "_rep[0-9]+\\.csv$")
  all_rep <- all_rep[grepl(base, basename(all_rep), ignore.case = FALSE)]
  if (!length(all_rep)) {
    log_line(sprintf("No replicate files matched --species %s; exiting.", opts$species), logf)
    h2o.shutdown(prompt = FALSE)
    quit(status = 0)
  }
}

# --- AutoML settings (same pattern as 03) ---
auto_seed <- cfg$h2o$seed %||% 123L
if (toupper(cfg$mode) == "REPRO") {
  include_algos    <- cfg$h2o$automl$include_algos_repro %||% c("GBM","DRF")
  max_models       <- opts$max_models %||% cfg$h2o$automl$max_models_repro %||% 20L
} else {
  include_algos    <- cfg$h2o$automl$include_algos_fast
  if (is.null(include_algos) || length(include_algos) == 0L) {
    include_algos  <- c("GBM","DRF","XGBoost")
  }
  max_models       <- opts$max_models %||% cfg$h2o$automl$max_models_fast %||% 60L
}
max_runtime_secs   <- opts$max_runtime_secs %||% 120L

# --- helpers ---
make_spatial_folds <- function(df, block_km = 2.5, k = 5L) {
  rngx <- range(df$lon); rngy <- range(df$lat)
  step <- block_km / 111
  bid_x <- floor((df$lon - rngx[1]) / max(step, .Machine$double.eps))
  bid_y <- floor((df$lat - rngy[1]) / max(step, .Machine$double.eps))
  block_id <- bid_x * 1e6 + bid_y
  set.seed( (cfg$seeds$cv_base %||% 22000L) )
  u <- sort(unique(block_id))
  fold_map <- setNames( (seq_along(u) - 1L) %% k + 1L, u)
  df$fold_id <- unname(fold_map[as.character(block_id)])
  df
}

drop_nzv <- function(df, cols) {
  keep <- vapply(cols, function(cl) {
    v <- df[[cl]]
    if (is.factor(v)) return(length(unique(v)) > 1L)
    suppressWarnings(sd(v, na.rm = TRUE)) > 0
  }, logical(1))
  cols[keep]
}

# deterministic species/rep seeds
seed_for_rep <- function(sp, rep_idx, base) {
  h <- sum(utf8ToInt(as.character(sp))) + as.integer(rep_idx) * 131
  as.integer(abs(h + base) %% .Machine$integer.max)
}

# --- process each replicate file ---
for (f in sort(all_rep)) {
  fname   <- basename(f)                 # e.g., "E5A_rep1.csv"
  sp      <- sub("_rep[0-9]+\\.csv$", "", fname)
  rep_idx <- as.integer(sub("^.*_rep([0-9]+)\\.csv$", "\\1", fname))

  out_dir <- file.path(paths$results_h2o, opts$run, sp, "replicates", paste0("rep", rep_idx))
  dir_ensure(out_dir)

  # resume guard
  if (skip_if_done(out_dir, c("leader_model.txt", "metrics_in_sample.csv", paste0("prediction_", sp, "_rep", rep_idx, ".tif")))) {
    log_line(sprintf("Skip %s (rep %d) — outputs present", sp, rep_idx), logf)
    next
  }

  log_line(sprintf("Training %s (rep %d)", sp, rep_idx), logf)

  # --- load replicate occurrence (presence only) ---
  sp_df <- readr::read_csv(f, show_col_types = FALSE)
  stopifnot(all(c("lon","lat") %in% names(sp_df)))
  pres <- sp_df %>% dplyr::select(lon, lat) %>% dplyr::mutate(pa = 1L)

  # --- deterministic background (depends on species + rep) ---
  bg_base <- (cfg$seeds$bg_base %||% 21000L)
  set.seed(seed_for_rep(sp, rep_idx, bg_base))
  bg_pts <- terra::spatSample(Renv[[1]], size = nrow(pres), method = "random", as.points = TRUE, na.rm = TRUE)
  bg_xy  <- terra::crds(bg_pts)
  bg_df  <- data.frame(lon = bg_xy[,1], lat = bg_xy[,2], pa = 0L)
  df_sp  <- dplyr::bind_rows(pres, bg_df)

  # --- extract predictors & prune ---
  pts  <- terra::vect(df_sp, geom = c("lon","lat"), crs = terra::crs(Renv))
  vals <- terra::extract(Renv, pts)[,-1, drop = FALSE]
  df   <- dplyr::bind_cols(df_sp, as.data.frame(vals)) %>% tidyr::drop_na()

  x_vars <- intersect(keep_vars, names(df))
  x_vars <- drop_nzv(df, x_vars)

  # guards
  n_min  <- 50L
  cls_tab <- table(df$pa)
  if (nrow(df) < n_min) {
    log_line(sprintf("Skip %s (rep %d) — too few rows after drop_na(): %d < %d", sp, rep_idx, nrow(df), n_min), logf); next
  }
  if (length(cls_tab) < 2 || any(cls_tab < 10)) {
    log_line(sprintf("Skip %s (rep %d) — class imbalance after drop_na() (counts: %s)",
                     sp, rep_idx, paste(names(cls_tab), as.integer(cls_tab), sep="=", collapse=", ")), logf); next
  }
  if (length(x_vars) == 0L) {
    log_line(sprintf("Skip %s (rep %d) — no usable predictors after pruning/NZV", sp, rep_idx), logf); next
  }

  # spatial folds
  df <- make_spatial_folds(df,
                           block_km = cfg$spatial_cv$block_km %||% 2.5,
                           k        = cfg$spatial_cv$folds     %||% 5L)

  # H2O frames
  hf <- as.h2o(df); hf["pa"] <- as.factor(hf[["pa"]])

  # --- AutoML with fallback ---
  aml <- tryCatch(
    h2o.automl(
      x = x_vars, y = "pa",
      training_frame   = hf,
      seed             = auto_seed + rep_idx,     # slight offset per replicate, still deterministic
      include_algos    = include_algos,
      max_models       = max_models,
      max_runtime_secs = max_runtime_secs,
      sort_metric      = "AUC",
      balance_classes  = TRUE
    ),
    error = function(e) {
      log_line(sprintf("AutoML error for %s (rep %d): %s", sp, rep_idx, conditionMessage(e)), logf)
      NULL
    }
  )

  leader <- NULL
  if (!is.null(aml)) {
    lb <- tryCatch(as.data.frame(aml@leaderboard), error = function(e) NULL)
    if (!is.null(lb) && nrow(lb) > 0) {
      leader <- aml@leader
    } else {
      log_line(sprintf("Empty AutoML leaderboard for %s (rep %d) — fallback GBM.", sp, rep_idx), logf)
    }
  } else {
    log_line(sprintf("AutoML returned NULL for %s (rep %d) — fallback GBM.", sp, rep_idx), logf)
  }

  if (is.null(leader)) {
    leader <- h2o.gbm(
      x = x_vars, y = "pa",
      training_frame = hf,
      distribution   = "bernoulli",
      ntrees         = 100,
      max_depth      = 5,
      learn_rate     = 0.05,
      seed           = auto_seed + rep_idx,
      balance_classes= TRUE
    )
  }

  # save model + info
  model_dir <- h2o.saveModel(leader, path = out_dir, force = TRUE)
  writeLines(c(
    paste0("leader_id=", leader@model_id),
    paste0("n_vars=", length(x_vars)),
    paste0("nthreads=", threads),
    paste0("mode=", mode),
    paste0("run=", opts$run),
    paste0("replicate=", rep_idx)
  ), file.path(out_dir, "leader_model.txt"))

  # in-sample metrics
  perf <- h2o.performance(leader, newdata = hf)
  ins  <- data.frame(
    dataset = sp, replicate = rep_idx,
    rmse    = tryCatch(h2o.rmse(perf),    error = function(e) NA_real_),
    auc     = tryCatch(h2o.auc(perf),     error = function(e) NA_real_),
    logloss = tryCatch(h2o.logloss(perf), error = function(e) NA_real_)
  )
  readr::write_csv(ins, file.path(out_dir, "metrics_in_sample.csv"))

  # varimp (if available)
  varimp <- tryCatch(as.data.frame(h2o.varimp(leader)), error = function(e) NULL)
  if (!is.null(varimp) && nrow(varimp) > 0) {
    readr::write_csv(varimp, file.path(out_dir, sprintf("varimp_%s_rep%d.csv", sp, rep_idx)))
    suppressWarnings({
      pdf(file.path(out_dir, sprintf("varimp_%s_rep%d.pdf", sp, rep_idx)))
      print(h2o.varimp_plot(leader, num_of_features = 10))
      dev.off()
    })
  }

  # partial dependence (first 5 vars)
  pp_vars <- head(x_vars, 5)
  suppressWarnings({
    pdf(file.path(out_dir, sprintf("partial_%s_rep%d.pdf", sp, rep_idx)))
    print(h2o.partialPlot(object = leader, data = hf, cols = pp_vars))
    dev.off()
  })

  # raster prediction (chunked)
  pred_r <- predict_raster_h2o(Renv[[x_vars]], leader, block_rows = 50000)
  pred_path <- file.path(out_dir, sprintf("prediction_%s_rep%d.tif", sp, rep_idx))
  terra::writeRaster(pred_r, filename = pred_path, overwrite = TRUE)

  log_line(sprintf("Done %s (rep %d)", sp, rep_idx), logf)
}

log_line("All replicate jobs finished.", logf)
h2o.shutdown(prompt = FALSE)
