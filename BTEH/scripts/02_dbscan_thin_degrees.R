#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(fs)
  library(dbscan)
  library(yaml)
  library(tibble)
})

# --- Locate script dir and project root, so sourcing works from any CWD ---
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

# ---------------- CLI ----------------
opt <- list(
  make_option(c("--run"),         type="character", default="B",     help="Run tag (e.g., A or B)"),
  make_option(c("--eps_deg"),     type="double",    default=0.004,   help="DBSCAN eps in *degrees* (~400 m at equator)"),
  make_option(c("--minPts"),      type="integer",   default=15,      help="DBSCAN minPts"),
  make_option(c("--fraction"),    type="double",    default=0.09,    help="Sampling fraction within clusters/noise"),
  make_option(c("--min_samples"), type="integer",   default=4,       help="Minimum samples to take per cluster"),
  make_option(c("--input"),       type="character", default=NULL,    help="Input folder of per-species cleaned CSVs (e.g., data/clean)"),
  make_option(c("--out_thin"),    type="character", default=NULL,    help="Output folder for MAIN replicate (default: data/occ/thinned_DBSCAN)"),
  make_option(c("--out_unc"),     type="character", default=NULL,    help="Output folder for UNCERTAINTY replicates (default: <project_root>/replicates_DBSCAN)"),
  make_option(c("--species"),     type="character", default=NULL,    help="Optional: process a single species (e.g., E3A)")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg  <- read_config()
mode <- set_mode(cfg)

# ---------------- Paths ----------------
in_dir    <- opts$input    %||% (cfg$paths$clean           %||% fs::path(.root, "data", "clean"))
out_thin  <- opts$out_thin %||% (cfg$paths$occ_thin        %||% fs::path(.root, "data","occ","thinned_DBSCAN"))
out_unc   <- opts$out_unc  %||% fs::path(.root, "replicates_DBSCAN")
logf      <- fs::path(cfg$paths$logs   %||% fs::path(.root, "logs"),  paste0("02_dbscan_thin_", opts$run, ".log"))
plans_yml <- fs::path(cfg$paths$plans  %||% fs::path(.root, "plans"), "dbscan_params.yml")
results_compare_dir <- cfg$paths$compare %||% fs::path(.root, "results", "compare")

dir_ensure(out_thin); dir_ensure(out_unc)
dir_ensure(fs::path_dir(logf)); dir_ensure(fs::path_dir(plans_yml))
dir_ensure(results_compare_dir)

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
files_df <- tibble(
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

# Return cleaned df + diagnostics for losses
normalize_with_diagnostics <- function(df, forced_species) {
  original <- df
  names(df) <- std_names(names(df))
  lon_candidates <- c("lon","longitude","x","long","long_wgs84","decimallongitude")
  lat_candidates <- c("lat","latitude","y","lat_wgs84","decimallatitude")
  pick <- function(cands) { hit <- intersect(cands, names(df)); if (length(hit)) hit[[1]] else NA_character_ }
  lon_col <- pick(lon_candidates); lat_col <- pick(lat_candidates)
  if (is.na(lon_col) || is.na(lat_col)) stop("Could not find longitude/latitude columns.")

  n_raw <- nrow(df)

  # Coerce to numeric and track parse NAs
  lon_num <- suppressWarnings(as.numeric(df[[lon_col]]))
  lat_num <- suppressWarnings(as.numeric(df[[lat_col]]))

  na_parse <- sum(!is.finite(lon_num) | !is.finite(lat_num))

  # Bounds filter
  in_bounds <- is.finite(lon_num) & is.finite(lat_num) &
               lon_num >= -180 & lon_num <= 180 &
               lat_num >= -90  & lat_num <= 90
  n_oob <- sum(!in_bounds & (is.finite(lon_num) & is.finite(lat_num)))

  # Keep only finite + in-bounds
  keep_idx <- which(in_bounds)
  kept_df <- tibble(species = forced_species, lon = lon_num[keep_idx], lat = lat_num[keep_idx])
  n_after_bounds <- nrow(kept_df)

  # Drop exact duplicate lon/lat
  before_distinct <- n_after_bounds
  kept_df <- kept_df %>% distinct(species, lon, lat, .keep_all = TRUE)
  n_dupes <- before_distinct - nrow(kept_df)

  out <- list(
    df = kept_df,
    stats = tibble(
      species = forced_species,
      n_raw = n_raw,
      n_na_or_nonfinite = na_parse,
      n_out_of_bounds = n_oob,
      n_after_clean = nrow(kept_df),
      n_exact_dupes_dropped = n_dupes
    )
  )
  return(out)
}

sample_with_dbscan <- function(df, eps_deg, minPts, fraction, min_samples, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  clu <- dbscan::dbscan(as.matrix(df[, c("lon","lat")]), eps = eps_deg, minPts = minPts)
  df$cluster_id <- as.integer(clu$cluster)
  df %>%
    group_by(cluster_id) %>%
    group_modify(~ {
      if (.y$cluster_id == 0L) {
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
total_raw <- 0L; total_clean <- 0L
total_main <- 0L; total_unc <- c(0L,0L,0L)
summaries <- list()
cleaning_stats <- list()

for (i in seq_len(nrow(files_df))) {
  sp   <- files_df$species[i]
  path <- files_df$path[i]

  raw <- safe_read_csv(path)

  # --- CLEAN & DIAGNOSE ---
  diag <- normalize_with_diagnostics(raw, sp)
  sp_df <- diag$df
  stat  <- diag$stats

  total_raw   <- total_raw   + stat$n_raw
  total_clean <- total_clean + stat$n_after_clean
  cleaning_stats[[length(cleaning_stats)+1L]] <- stat

  n_in_clean <- nrow(sp_df)

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
  log_line(sprintf("MAIN  %s: raw=%d | clean=%d -> kept=%d (saved %s)",
                   sp, stat$n_raw, n_in_clean, nrow(main_df), main_out_path), logf)

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
    total_unc[rep_idx]  <- total_unc[rep_idx] + nrow(rep_df)
    rep_out_path <- fs::path(out_unc, paste0(sp, "_rep", rep_idx, ".csv"))
    safe_write_csv(rep_df, rep_out_path)
    log_line(sprintf("UNC%1d %s: raw=%d | clean=%d -> kept=%d (saved %s)",
                     rep_idx, sp, stat$n_raw, n_in_clean, nrow(rep_df), rep_out_path), logf)
  }

  # per-species summary row
  summaries[[length(summaries)+1L]] <- tibble(
    species     = sp,
    input_raw   = stat$n_raw,
    input_clean = n_in_clean,
    kept_main   = nrow(main_df),
    kept_rep1   = rep_counts[1],
    kept_rep2   = rep_counts[2],
    kept_rep3   = rep_counts[3]
  )
}

# ------------- Wrap up: diagnostics + summaries -------------
clean_df <- dplyr::bind_rows(cleaning_stats)
sumdf    <- dplyr::bind_rows(summaries)

# Write detailed cleaning diagnostics per species
clean_diag_path <- fs::path(fs::path_dir(logf), paste0("cleaning_diagnostics_", opts$run, ".csv"))
safe_write_csv(clean_df, clean_diag_path)
log_line(paste0("Cleaning diagnostics saved: ", clean_diag_path), logf)

# Write thinning summary (before/after), for the front-end & QA
sumdf <- sumdf %>% mutate(
  keep_ratio_main = ifelse(input_clean > 0, kept_main / input_clean, NA_real_)
)
thin_summary_path <- fs::path(results_compare_dir, "thin_summary.csv")
safe_write_csv(sumdf, thin_summary_path)
log_line(paste0("Thinning summary saved: ", thin_summary_path), logf)

# Totals (raw vs clean vs kept)
totals <- tibble(
  run = opts$run,
  total_raw = total_raw,
  total_clean = total_clean,
  main_kept = total_main,
  unc_rep1_kept = total_unc[1],
  unc_rep2_kept = total_unc[2],
  unc_rep3_kept = total_unc[3],
  raw_minus_clean = total_raw - total_clean
)
totals_path <- fs::path(fs::path_dir(logf), paste0("thin_totals_", opts$run, ".csv"))
safe_write_csv(totals, totals_path)
log_line(paste0("Totals saved: ", totals_path), logf)

# Console print: neat table
cat("\n==================== DBSCAN THINNING SUMMARY ====================\n")
print(
  sumdf %>%
    select(species, input_raw, input_clean, kept_main, kept_rep1, kept_rep2, kept_rep3, keep_ratio_main) %>%
    arrange(species)
)
cat("-----------------------------------------------------------------\n")
print(totals)
cat("-----------------------------------------------------------------\n")
cat("Note: input_raw - input_clean = rows dropped during cleaning (NA/non-finite, out-of-bounds, or exact lon/lat duplicates).\n\n")

# Final log line
log_line(sprintf("DONE. Total RAW: %d | CLEAN: %d | main kept: %d | unc reps kept: [%d, %d, %d]",
                 total_raw, total_clean, total_main, total_unc[1], total_unc[2], total_unc[3]), logf)
