#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(fs)
  library(dbscan)
  library(yaml)
})

# --- Locate script dir and project root, so sourcing works from any CWD ---
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  filearg <- grep("^--file=", args, value = TRUE)
  if (length(filearg)) return(normalizePath(sub("^--file=", "", filearg)))
  # Fallback when run via source() or in IDE
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

# ---------------- CLI ----------------
opt <- list(
  make_option(c("--run"),         type="character", default="B",     help="Run tag (e.g., A or B)"),
  make_option(c("--eps_deg"),     type="double",    default=0.004,   help="DBSCAN eps in *degrees* (default 0.004 ≈ 400m at equator)"),
  make_option(c("--minPts"),      type="integer",   default=15,      help="DBSCAN minPts"),
  make_option(c("--fraction"),    type="double",    default=0.09,    help="Sampling fraction within clusters/noise"),
  make_option(c("--min_samples"), type="integer",   default=4,       help="Minimum samples to take per cluster"),
  make_option(c("--input"),       type="character", default=NULL,    help="Input folder of per-species cleaned CSVs (e.g., data/clean)"),
  make_option(c("--out_thin"),    type="character", default=NULL,    help="Output folder for MAIN replicate (default: data/occ/thinned_DBSCAN)"),
  # NOTE: default replicates folder is now at PROJECT ROOT for convenience (can be overridden)
  make_option(c("--out_unc"),     type="character", default=NULL,    help="Output folder for UNCERTAINTY replicates (default: <project_root>/replicates_DBSCAN)"),
  make_option(c("--species"),     type="character", default=NULL,    help="Optional: process a single species (e.g., E3A)")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg  <- read_config()
mode <- set_mode(cfg)

# ---------------- Paths ----------------
in_dir    <- opts$input    %||% (cfg$paths$clean           %||% fs::path(.root, "data", "clean"))
out_thin  <- opts$out_thin %||% (cfg$paths$occ_thin        %||% fs::path(.root, "data","occ","thinned_DBSCAN"))
# Default UNCERTAINTY/replicates folder at *project root* (overridable via --out_unc)
out_unc   <- opts$out_unc  %||% fs::path(.root, "replicates_DBSCAN")

logf      <- fs::path(cfg$paths$logs   %||% fs::path(.root, "logs"),  paste0("02_dbscan_thin_", opts$run, ".log"))
plans_yml <- fs::path(cfg$paths$plans  %||% fs::path(.root, "plans"), "dbscan_params.yml")

dir_ensure(out_thin)
dir_ensure(out_unc)
dir_ensure(fs::path_dir(logf))
dir_ensure(fs::path_dir(plans_yml))

log_line(paste0("Starting 02_dbscan_thin.R (mode=", mode, ")"), logf)
log_line(paste0("Run: ", opts$run,
                " | eps_deg: ", opts$eps_deg,
                " | minPts: ", opts$minPts,
                " | fraction: ", opts$fraction,
                " | min_samples: ", opts$min_samples), logf)
log_line(paste0("Input: ", in_dir,
                " | MAIN out: ", out_thin,
                " | UNCERTAINTY (replicates) out: ", out_unc), logf)

# ------------- Inputs -------------
csvs <- fs::dir_ls(in_dir, type = "file", recurse = FALSE)
csvs <- csvs[grepl("\\.csv$", tolower(csvs))]
if (!length(csvs)) {
  log_line("No CSV files found in input; exiting.", logf)
  quit(status = 0)
}

files_df <- tibble::tibble(
  path    = csvs,
  species = fs::path_ext_remove(fs::path_file(csvs))
)

if (!is.null(opts$species)) {
  files_df <- dplyr::filter(files_df, species == opts$species)
  if (!nrow(files_df)) {
    log_line(paste0("Species '", opts$species, "' not found in input; exiting."), logf)
    quit(status = 0)
  }
}

# ------------- Helpers -------------
std_names <- function(nms) {
  nms |>
    tolower() |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("\\.", "_")
}

normalize_columns_force_species <- function(df, forced_species) {
  names(df) <- std_names(names(df))
  lon_candidates <- c("lon","longitude","x","long","long_wgs84","decimallongitude")
  lat_candidates <- c("lat","latitude","y","lat_wgs84","decimallatitude")
  pick <- function(cands) { hit <- intersect(cands, names(df)); if (length(hit)) hit[[1]] else NA_character_ }
  lon_col <- pick(lon_candidates); lat_col <- pick(lat_candidates)
  if (is.na(lon_col) || is.na(lat_col)) stop("Could not find longitude/latitude columns.")
  df |>
    transmute(
      species = forced_species,
      lon = as.numeric(.data[[lon_col]]),
      lat = as.numeric(.data[[lat_col]])
    ) |>
    filter(is.finite(lon), is.finite(lat), lon >= -180, lon <= 180, lat >= -90, lat <= 90) |>
    distinct(species, lon, lat, .keep_all = TRUE)
}

# Selection logic (STRICTLY in *degrees*):
# - DBSCAN on (lon, lat) with eps = opts$eps_deg, minPts = opts$minPts
# - Noise (cluster==0): keep proportion = fraction, but at least 1 when there is >1 row (your original behavior)
# - Clusters (>0): keep max(min_samples, round(n * fraction)), capped at n
sample_with_dbscan <- function(df, eps_deg, minPts, fraction, min_samples, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  clu <- dbscan::dbscan(as.matrix(df[, c("lon","lat")]), eps = eps_deg, minPts = minPts)
  df$cluster_id <- as.integer(clu$cluster)
  df %>%
    group_by(cluster_id) %>%
    group_modify(~ {
      if (.y$cluster_id == 0L) {
        # noise: proportion (preserve your 'keep >=1' intent when >1 rows)
        if (nrow(.x) <= 1) return(.x)
        keep_n <- max(1, round(nrow(.x) * fraction))
        dplyr::slice_sample(.x, n = min(nrow(.x), keep_n))
      } else {
        n_to_sample <- max(min_samples, round(nrow(.x) * fraction))
        dplyr::slice_sample(.x, n = min(nrow(.x), n_to_sample))
      }
    }) %>%
    ungroup() %>%
    select(-cluster_id)
}

# ------------- Record params -------------
params <- list()
if (file.exists(plans_yml)) {
  params <- tryCatch(yaml::read_yaml(plans_yml), error = function(e) list())
}
if (is.null(params[[opts$run]])) params[[opts$run]] <- list()
params[[opts$run]]$dbscan <- list(
  eps_deg      = opts$eps_deg,
  minPts       = opts$minPts,
  fraction     = opts$fraction,
  min_samples  = opts$min_samples,
  input        = fs::path_abs(in_dir),
  main_out     = fs::path_abs(out_thin),
  uncertainty  = fs::path_abs(out_unc),
  timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
ok <- tryCatch({ yaml::write_yaml(params, plans_yml); TRUE }, error = function(e) { FALSE })
if (!ok) log_line(paste0("WARN: failed to write params YAML: ", plans_yml), logf) else
  log_line(paste0("Recorded params in: ", plans_yml), logf)

# ------------- Process -------------
total_in <- 0L; total_main <- 0L; total_unc <- c(0L,0L,0L)
summaries <- list()

for (i in seq_len(nrow(files_df))) {
  sp   <- files_df$species[i]
  path <- files_df$path[i]

  raw <- safe_read_csv(path)
  sp_df <- normalize_columns_force_species(raw, sp)
  n_in  <- nrow(sp_df)
  total_in <- total_in + n_in

  # Seeds: deterministic for main; fixed replicates 1..3 but derived per species
  main_seed <- seed_for(paste0("dbscan_main_", opts$run, "_", sp),
                        base = cfg$seeds$per_dataset_base %||% 30000L)

  # ---- MAIN replicate (deterministic) ----
  main_df <- sample_with_dbscan(
    df = sp_df,
    eps_deg = opts$eps_deg,
    minPts = opts$minPts,
    fraction = opts$fraction,
    min_samples = opts$min_samples,
    seed = main_seed
  )
  total_main <- total_main + nrow(main_df)
  main_out_path <- fs::path(out_thin, paste0(sp, ".csv"))
  safe_write_csv(main_df, main_out_path)
  log_line(sprintf("MAIN  %s: %d -> %d (saved %s)", sp, n_in, nrow(main_df), main_out_path), logf)

  # ---- UNCERTAINTY replicates (rep1..rep3) ----
  rep_counts <- integer(3)
  for (rep_idx in 1:3) {
    rep_seed <- seed_for(paste0("dbscan_unc_", opts$run, "_", sp, "_rep", rep_idx),
                         base = (cfg$seeds$per_dataset_base %||% 30000L) + rep_idx)
    rep_df <- sample_with_dbscan(
      df = sp_df,
      eps_deg = opts$eps_deg,
      minPts = opts$minPts,
      fraction = opts$fraction,
      min_samples = opts$min_samples,
      seed = rep_seed
    )
    rep_counts[rep_idx] <- nrow(rep_df)
    total_unc[rep_idx] <- total_unc[rep_idx] + nrow(rep_df)
    rep_out_path <- fs::path(out_unc, paste0(sp, "_rep", rep_idx, ".csv"))
    safe_write_csv(rep_df, rep_out_path)
    log_line(sprintf("UNC%1d %s: %d -> %d (saved %s)", rep_idx, sp, n_in, nrow(rep_df), rep_out_path), logf)
  }

  # optional per-species summary row
  summaries[[length(summaries)+1L]] <- tibble::tibble(
    species   = sp,
    n_input   = n_in,
    kept_main = nrow(main_df),
    kept_rep1 = rep_counts[1],
    kept_rep2 = rep_counts[2],
    kept_rep3 = rep_counts[3]
  )
}

# ------------- Wrap up -------------
if (length(summaries)) {
  sumdf <- dplyr::bind_rows(summaries)
  sum_path <- fs::path(fs::path_dir(logf), paste0("dbscan_summary_", opts$run, ".csv"))
  safe_write_csv(sumdf, sum_path)
  log_line(paste0("Summary saved: ", sum_path), logf)
}

log_line(sprintf("DONE. Total input: %d | main kept: %d | unc reps kept: [%d, %d, %d]",
                 total_in, total_main, total_unc[1], total_unc[2], total_unc[3]), logf)
