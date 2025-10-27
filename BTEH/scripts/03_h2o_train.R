# ============================================================================
# Purpose: Optimized training pipeline — modular, resumable, deterministic in
#          REPRO mode; FAST mode uses all cores minus 2. Reads thinned points
#          from data/occ/thinned_DBSCAN/*.csv and saves to results/H2O/<RUN>/<SP>/.
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
  make_option(c("--run"),              type = "character", default = "B",   help = "Run tag: A or B"),
  make_option(c("--mode"),             type = "character", default = NULL,  help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"),          type = "character", default = NULL,  help = "Optional single species (e.g., E3A)"),
  make_option(c("--max_models"),       type = "integer",   default = NULL,  help = "Override AutoML max_models"),
  make_option(c("--max_runtime_secs"), type = "integer",   default = NULL,  help = "Optional AutoML time budget (seconds)")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)

# --- paths ---
occ_dir   <- cfg$paths$occ %||% file.path("data", "occ", "thinned_DBSCAN")
res_root  <- cfg$paths$results_h2o %||% file.path("results", "H2O")
plans_dir <- cfg$paths$plans       %||% "plans"
logs_dir  <- cfg$paths$logs        %||% "logs"
logf      <- file.path(logs_dir, sprintf("03_h2o_train_%s.log", opts$run))

dir_ensure(res_root); dir_ensure(plans_dir); dir_ensure(logs_dir)
log_line(sprintf("Starting 03_h2o_train.R  (mode=%s, run=%s)", mode, opts$run), logf)

# --- choose raster dir by run tag ---
if (toupper(opts$run) == "A") {
  env_dir <- cfg$paths$envi_after %||% file.path("data","envi","A")
} else {
  env_dir <- cfg$paths$envi_before %||% file.path("data","envi","B")
}
ras_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(ras_files) > 0)
ras_names <- make.names(tools::file_path_sans_ext(basename(ras_files)), unique = TRUE)
Renv <- terra::rast(ras_files); names(Renv) <- ras_names

# --- Kendall plan (read or create once per run) ---
keepvars_csv <- file.path(plans_dir, opts$run, "keepvars.csv")
if (file.exists(keepvars_csv)) {
  keep_vars <- readr::read_csv(keepvars_csv, show_col_types = FALSE)$variable
  log_line(sprintf("Read keepvars: %d variables", length(keep_vars)), logf)
} else {
  keep_vars <- kendall_plan(opts$run, Renv,
                            out_dir = file.path(plans_dir, opts$run),
                            cutoff   = cfg$kendall$cutoff %||% 0.8,
                            nsample  = cfg$kendall$sample %||% 5000L)
  log_line(sprintf("Created keepvars: %d variables", length(keep_vars)), logf)
}
Renv <- Renv[[keep_vars]]

# --- H2O init with proper threads ---
threads <- init_h2o(cfg)  # utils_h2o.R: no h2o.setSeed (OK for H2O 3.46+)
log_line(sprintf("H2O threads: %d", threads), logf)

# --- species list ---
occ_files <- list.files(occ_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(occ_files) > 0)
sp_list <- sort(tools::file_path_sans_ext(basename(occ_files)))
if (!is.null(opts$species)) sp_list <- intersect(sp_list, opts$species)

# --- AutoML settings from config ---
auto_seed <- cfg$h2o$seed %||% 123L
if (toupper(cfg$mode) == "REPRO") {
  include_algos <- cfg$h2o$automl$include_algos_repro %||% c("GBM","DRF")
  max_models    <- opts$max_models %||% cfg$h2o$automl$max_models_repro %||% 20L
} else {
  include_algos <- cfg$h2o$automl$include_algos_fast
  if (is.null(include_algos) || length(include_algos) == 0L) {
    include_algos <- c("GBM","DRF","XGBoost")  # ensure at least a few fast learners
  }
  max_models    <- opts$max_models %||% cfg$h2o$automl$max_models_fast %||% 60L
}
max_runtime_secs <- opts$max_runtime_secs %||% 120L  # default small budget for smoke tests

# --- deterministic fold generator (spatial blocks) ---
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

# --- helper: drop zero-variance predictors ---
drop_nzv <- function(df, cols) {
  keep <- vapply(cols, function(cl) {
    v <- df[[cl]]
    if (is.factor(v)) return(length(unique(v)) > 1L)
    suppressWarnings(sd(v, na.rm = TRUE)) > 0
  }, logical(1))
  cols[keep]
}

for (sp in sp_list) {
  sp_csv <- file.path(occ_dir, paste0(sp, ".csv"))
  if (!file.exists(sp_csv)) next

  out_dir <- file.path(res_root, opts$run, sp)
  dir_ensure(out_dir)

  # Resume guard
  if (skip_if_done(out_dir, c("leader_model.txt", "metrics_in_sample.csv"))) {
    log_line(sprintf("Skip %s — outputs present", sp), logf)
    next
  }

  log_line(sprintf("Training %s", sp), logf)

  sp_df <- readr::read_csv(sp_csv, show_col_types = FALSE)
  stopifnot(all(c("lon","lat") %in% names(sp_df)))

  # Presence/absence with deterministic background
  pres <- sp_df %>% dplyr::select(lon, lat) %>% dplyr::mutate(pa = 1L)
  sp_hash <- sum(utf8ToInt(as.character(sp))) %% 10000L
  set.seed( (cfg$seeds$bg_base %||% 21000L) + sp_hash )
  bg_pts <- terra::spatSample(Renv[[1]], size = nrow(pres), method = "random",
                              as.points = TRUE, na.rm = TRUE)
  bg_xy  <- terra::crds(bg_pts)
  bg_df  <- data.frame(lon = bg_xy[,1], lat = bg_xy[,2], pa = 0L)
  df_sp  <- dplyr::bind_rows(pres, bg_df)

  # Extract predictors (drop NAs)
  pts  <- terra::vect(df_sp, geom = c("lon","lat"), crs = terra::crs(Renv))
  vals <- terra::extract(Renv, pts)[,-1, drop = FALSE]
  df   <- dplyr::bind_cols(df_sp, as.data.frame(vals)) %>% tidyr::drop_na()

  # Predictors intersection + NZV filter
  x_vars <- intersect(keep_vars, names(df))
  x_vars <- drop_nzv(df, x_vars)

  # Guardrails
  n_min  <- 50L
  cls_tab <- table(df$pa)
  if (nrow(df) < n_min) {
    log_line(sprintf("Skip %s — too few rows after drop_na(): %d < %d", sp, nrow(df), n_min), logf)
    next
  }
  if (length(cls_tab) < 2 || any(cls_tab < 10)) {
    log_line(sprintf("Skip %s — class imbalance after drop_na() (counts: %s)", sp,
                     paste(names(cls_tab), as.integer(cls_tab), sep="=", collapse=", ")), logf)
    next
  }
  if (length(x_vars) == 0L) {
    log_line(sprintf("Skip %s — no usable predictors after pruning/NZV", sp), logf)
    next
  }

  # Spatial folds (deterministic grid blocks)
  df <- make_spatial_folds(df,
                           block_km = cfg$spatial_cv$block_km %||% 2.5,
                           k        = cfg$spatial_cv$folds     %||% 5L)

  # H2O frames
  hf <- as.h2o(df); hf["pa"] <- as.factor(hf[["pa"]])

  # --- AutoML with safe fallback ---
  aml <- tryCatch(
    h2o.automl(
      x = x_vars, y = "pa",
      training_frame   = hf,
      seed             = auto_seed,
      include_algos    = include_algos,
      max_models       = max_models,
      max_runtime_secs = max_runtime_secs,
      sort_metric      = "AUC",
      balance_classes  = TRUE
    ),
    error = function(e) {
      log_line(sprintf("AutoML error for %s: %s", sp, conditionMessage(e)), logf)
      NULL
    }
  )

  leader <- NULL
  if (!is.null(aml)) {
    # check leaderboard is non-empty
    lb <- tryCatch(as.data.frame(aml@leaderboard), error = function(e) NULL)
    if (!is.null(lb) && nrow(lb) > 0) {
      leader <- aml@leader
    } else {
      log_line(sprintf("Empty AutoML leaderboard for %s — switching to fallback GBM.", sp), logf)
    }
  } else {
    log_line(sprintf("AutoML returned NULL for %s — switching to fallback GBM.", sp), logf)
  }

  # Fallback: simple GBM (guarantees a model)
  if (is.null(leader)) {
    leader <- h2o.gbm(
      x = x_vars, y = "pa",
      training_frame = hf,
      distribution   = "bernoulli",
      ntrees         = 100,
      max_depth      = 5,
      learn_rate     = 0.05,
      seed           = auto_seed,
      balance_classes= TRUE
    )
  }

  # Save leader model + info
  model_dir <- h2o.saveModel(leader, path = out_dir, force = TRUE)
  writeLines(c(
    paste0("leader_id=", leader@model_id),
    paste0("n_vars=", length(x_vars)),
    paste0("nthreads=", threads),
    paste0("mode=", mode)
  ), file.path(out_dir, "leader_model.txt"))

  # In-sample performance (reference)
  perf <- h2o.performance(leader, newdata = hf)
  ins  <- data.frame(
    dataset = sp,
    rmse    = tryCatch(h2o.rmse(perf),    error = function(e) NA_real_),
    auc     = tryCatch(h2o.auc(perf),     error = function(e) NA_real_),
    logloss = tryCatch(h2o.logloss(perf), error = function(e) NA_real_)
  )
  readr::write_csv(ins, file.path(out_dir, "metrics_in_sample.csv"))

  # Variable importance (if available)
  varimp <- tryCatch(as.data.frame(h2o.varimp(leader)), error = function(e) NULL)
  if (!is.null(varimp) && nrow(varimp) > 0) {
    readr::write_csv(varimp, file.path(out_dir, sprintf("varimp_%s.csv", sp)))
    suppressWarnings({
      pdf(file.path(out_dir, sprintf("varimp_%s.pdf", sp)));
      print(h2o.varimp_plot(leader, num_of_features = 10));
      dev.off()
    })
  }

  # Partial dependence for first few variables
  pp_vars <- head(x_vars, 5)
  suppressWarnings({
    pdf(file.path(out_dir, sprintf("partial_%s.pdf", sp)));
    print(h2o.partialPlot(object = leader, data = hf, cols = pp_vars));
    dev.off()
  })

  # Raster prediction (chunked)
  pred_r <- predict_raster_h2o(Renv[[x_vars]], leader, block_rows = 50000)
  terra::writeRaster(pred_r,
                     filename = file.path(out_dir, sprintf("prediction_%s.tif", sp)),
                     overwrite = TRUE)

  log_line(sprintf("Done %s", sp), logf)
}

log_line("All species finished.", logf)
h2o.shutdown(prompt = FALSE)
