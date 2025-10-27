# =============================================================================
# 05_h2o_vs_ssdm_results.R — Between-method comparison from *replicates only*
# Sources:
#   - H2O:  results/H2O/<RUN>/<SP>/replicates/rep<k>/pred.tif
#   - SSDM: results/SSDM/<RUN>/<SP>/replicates/rep<k>/ESDM_<SP>_rep<k>.tif
# Outputs:
#   - results/compare/<RUN>/01_between_methods/{rasters,plots}/...
#   - results/compare/<RUN>/02_tables/{per_rep_metrics.csv, per_species_summary.csv}
#   - results/compare/<RUN>/03_maps/{<SP>_H2O_mean.png, <SP>_SSDM_mean.png, diff_mean.png}
# Notes:
#   - REPRO/FAST respected via config.yml
#   - NA-safe; resumes fine if files exist
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tools)
  library(grid)    # unit() in legends
  library(purrr)
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

source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- CLI ---------------------------------------------------------------------
opt <- list(
  make_option(c("--run"),     type = "character", default = "B",   help = "Run tag: A or B"),
  make_option(c("--mode"),    type = "character", default = NULL,  help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"), type = "character", default = NULL,  help = "Optional single species (e.g., E5A)"),
  make_option(c("--q"),       type = "double",    default = 0.75,  help = "Hotspot quantile (default 0.75)")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)

# --- paths & logs ------------------------------------------------------------
h2o_root     <- cfg$paths$results_h2o %||% file.path("results","H2O")
ssdm_root    <- cfg$paths$results_ssdm %||% file.path("results","SSDM")
compare_root <- file.path("results","compare", toupper(opts$run))
logs_dir     <- cfg$paths$logs %||% "logs"
dir_ensure(compare_root); dir_ensure(logs_dir)
dir_ensure(file.path(compare_root, "01_between_methods/rasters"))
dir_ensure(file.path(compare_root, "01_between_methods/plots"))
dir_ensure(file.path(compare_root, "02_tables"))
dir_ensure(file.path(compare_root, "03_maps"))

logf <- file.path(logs_dir, sprintf("05_compare_replicates_%s.log", toupper(opts$run)))
log_line(sprintf("Starting 05_h2o_vs_ssdm_results.R (mode=%s, run=%s, q=%.2f)",
                 cfg$mode, toupper(opts$run), opts$q), logf)

# --- helpers -----------------------------------------------------------------
align_to <- function(r1, r2, categorical = FALSE) {
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- project(r2, crs(r1), method = if (categorical) "near" else "bilinear")
    r2 <- resample(r2, r1,      method = if (categorical) "near" else "bilinear")
    r2 <- crop(r2, r1)
  }
  m <- !is.na(r1) & !is.na(r2)
  r1 <- mask(r1, m, maskvalues = 0)
  r2 <- mask(r2, m, maskvalues = 0)
  list(r1 = r1, r2 = r2)
}

jaccard_binary <- function(b1, b2) {
  inter <- global(b1 & b2, "sum", na.rm = TRUE)[[1]]
  union <- global(b1 | b2, "sum", na.rm = TRUE)[[1]]
  ifelse(union == 0, NA_real_, inter / union)
}

plot_raster_continuous <- function(
  r, title, out_png,
  center0 = FALSE,
  legend_title = "Suitability",
  barheight_cm = 12,
  barwidth_cm  = 0.7,
  base_size    = 12
) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() + coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = legend_title) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text  = element_text(size = base_size)
    )
  if (center0) {
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p <- p + scale_fill_gradientn(
      colors = c("#2c7fb8","#ffffbf","#d7191c"),
      limits = c(-lim, lim),
      na.value = NA,
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth  = unit(barwidth_cm,  "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  } else {
    p <- p + scale_fill_viridis_c(
      na.value = NA,
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth  = unit(barwidth_cm,  "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  }
  ggsave(out_png, p, width = 8.2, height = 6.2, dpi = 300)
}

# Discover species + replicate indices present in BOTH methods -----------------
run_tag <- toupper(opts$run)
h2o_run_dir  <- file.path(h2o_root,  run_tag)
ssdm_run_dir <- file.path(ssdm_root, run_tag)

species_h2o  <- list.dirs(h2o_run_dir,  full.names = FALSE, recursive = FALSE)
species_ssdm <- list.dirs(ssdm_run_dir, full.names = FALSE, recursive = FALSE)

if (!is.null(opts$species)) {
  species_h2o  <- intersect(species_h2o,  opts$species)
  species_ssdm <- intersect(species_ssdm, opts$species)
}
species_both <- intersect(species_h2o, species_ssdm)
stopifnot(length(species_both) > 0)
log_line(sprintf("Species considered (replicate outputs only): %s",
                 paste(species_both, collapse = ", ")), logf)

# Functions to list replicate rasters
list_h2o_rep_preds <- function(sp) {
  rep_dir <- file.path(h2o_run_dir, sp, "replicates")
  reps <- list.dirs(rep_dir, full.names = FALSE, recursive = FALSE)
  if (!length(reps)) return(tibble(rep = integer(), path = character()))
  tibble(
    rep  = as.integer(sub("^rep", "", reps)),
    path = file.path(rep_dir, reps, "pred.tif")
  ) %>% filter(file.exists(path))
}

list_ssdm_rep_preds <- function(sp) {
  rep_dir <- file.path(ssdm_run_dir, sp, "replicates")
  reps <- list.dirs(rep_dir, full.names = FALSE, recursive = FALSE)
  if (!length(reps)) return(tibble(rep = integer(), path = character()))
  tibble(
    rep  = as.integer(sub("^rep", "", reps)),
    path = file.path(rep_dir, reps, sprintf("ESDM_%s_rep%d.tif", sp, as.integer(sub("^rep","",reps))))
  ) %>% filter(file.exists(path))
}

# Main compute loop -----------------------------------------------------------
per_rep_rows <- list()

for (sp in species_both) {
  h2o_tbl  <- list_h2o_rep_preds(sp)
  ssdm_tbl <- list_ssdm_rep_preds(sp)

  reps_both <- intersect(h2o_tbl$rep, ssdm_tbl$rep)
  if (!length(reps_both)) {
    log_line(sprintf("Skip %s — no replicate overlap between H2O and SSDM.", sp), logf)
    next
  }

  log_line(sprintf("Processing %s | shared replicates: %s",
                   sp, paste(sort(reps_both), collapse = ", ")), logf)

  # Containers to compute mean rasters over shared reps
  diffs  <- list()
  predsH <- list()
  predsS <- list()

  for (rk in sort(reps_both)) {
    f_h2o  <- h2o_tbl$path[h2o_tbl$rep == rk][1]
    f_ssdm <- ssdm_tbl$path[ssdm_tbl$rep == rk][1]

    r_h2o  <- rast(f_h2o)
    r_ssdm <- rast(f_ssdm)

    al     <- align_to(r_ssdm, r_h2o, categorical = FALSE)
    r_ssdm <- al$r1; r_h2o <- al$r2

    v1 <- values(r_ssdm, mat = FALSE)
    v2 <- values(r_h2o,  mat = FALSE)
    idx <- !is.na(v1) & !is.na(v2)

    pearson  <- suppressWarnings(cor(v1[idx], v2[idx], method = "pearson"))
    spearman <- suppressWarnings(cor(v1[idx], v2[idx], method = "spearman"))
    rmse     <- sqrt(mean((v1[idx] - v2[idx])^2))
    mae      <- mean(abs(v1[idx] - v2[idx]))

    # Hotspot Jaccard at quantile q per replicate
    q1 <- quantile(v1[idx], opts$q, na.rm = TRUE)
    q2 <- quantile(v2[idx], opts$q, na.rm = TRUE)
    h1 <- r_ssdm > q1
    h2 <- r_h2o  > q2
    jac_hot <- jaccard_binary(h1, h2)

    per_rep_rows[[length(per_rep_rows) + 1L]] <- tibble(
      run = run_tag, dataset = sp, replicate = rk,
      pearson = pearson, spearman = spearman,
      rmse = rmse, mae = mae, jaccard_hotspot_q = jac_hot
    )

    # Collect for means
    predsH[[length(predsH) + 1L]] <- r_h2o
    predsS[[length(predsS) + 1L]] <- r_ssdm
    diffs[[length(diffs) + 1L]]   <- (r_h2o - r_ssdm)
  }

  # Write mean maps for this species (over shared reps)
  rH_mean <- do.call(c, predsH) |> mean(na.rm = TRUE)
  rS_mean <- do.call(c, predsS) |> mean(na.rm = TRUE)
  d_mean  <- do.call(c, diffs)  |> mean(na.rm = TRUE)

  writeRaster(rH_mean,
              file.path(compare_root, "01_between_methods/rasters",
                        sprintf("%s_H2O_mean.tif", sp)),
              overwrite = TRUE, datatype = "FLT4S",
              gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
  writeRaster(rS_mean,
              file.path(compare_root, "01_between_methods/rasters",
                        sprintf("%s_SSDM_mean.tif", sp)),
              overwrite = TRUE, datatype = "FLT4S",
              gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
  writeRaster(d_mean,
              file.path(compare_root, "01_between_methods/rasters",
                        sprintf("diff_mean_H2O_minus_SSDM_%s.tif", sp)),
              overwrite = TRUE, datatype = "FLT4S",
              gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))

  # Plots with bigger/longer legends
  plot_raster_continuous(rS_mean, paste0(sp, " — SSDM mean suitability"),
                         file.path(compare_root, "03_maps", sprintf("%s_SSDM_mean.png", sp)))
  plot_raster_continuous(rH_mean, paste0(sp, " — H2O mean suitability"),
                         file.path(compare_root, "03_maps", sprintf("%s_H2O_mean.png", sp)))
  plot_raster_continuous(d_mean,  paste0("Between-method mean difference (H2O − SSDM): ", sp),
                         file.path(compare_root, "01_between_methods/plots",
                                   sprintf("diff_mean_H2O_minus_SSDM_%s.png", sp)),
                         center0 = TRUE, legend_title = "Δ suitability")
}

# Per-replicate table ----------------------------------------------------------
per_rep_df <- bind_rows(per_rep_rows)
out_per_rep <- file.path(compare_root, "02_tables", "per_rep_metrics.csv")
readr::write_csv(per_rep_df, out_per_rep)
log_line(sprintf("Wrote per-replicate metrics: %s", out_per_rep), logf)

# Per-species summary (means, sds, n) -----------------------------------------
summary_df <- per_rep_df %>%
  group_by(run, dataset) %>%
  summarise(
    n_reps = dplyr::n(),
    across(c(pearson, spearman, rmse, mae, jaccard_hotspot_q),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  ) %>% ungroup()

out_summary <- file.path(compare_root, "02_tables", "per_species_summary.csv")
readr::write_csv(summary_df, out_summary)
log_line(sprintf("Wrote per-species summary: %s", out_summary), logf)

# --- annotated bar charts from per-species MEANS ------------------------------
bm_df <- summary_df %>%
  transmute(
    dataset,
    pearson        = pearson_mean,
    spearman       = spearman_mean,
    rmse           = rmse_mean,
    mae            = mae_mean,
    jaccard_hotspot_q = jaccard_hotspot_q_mean
  )

mpl_blue <- "#1f77b4"
# Order that looks nice; only keep present species
target_order   <- c("E1B","E5B","E5A","E2B","E4A","E3B","E3A","E4B","E6B")
present_order  <- intersect(target_order, bm_df$dataset)
if (length(present_order) == 0) present_order <- sort(unique(bm_df$dataset))
bm_df <- bm_df %>% mutate(dataset = factor(dataset, levels = rev(present_order))) %>% arrange(dataset)

save_bar <- function(df, metric, title, fname, clamp01 = TRUE, digits = 3) {
  df <- df %>% mutate(label = sprintf(paste0("%.", digits, "f"), .data[[metric]]))
  maxv <- suppressWarnings(max(df[[metric]], na.rm = TRUE))
  upper <- if (clamp01) 1 else max(1e-9, maxv) * 1.10

  p <- ggplot(df, aes(x = dataset, y = .data[[metric]])) +
    geom_col(fill = mpl_blue, width = 0.8) +
    geom_text(aes(label = label), hjust = -0.15, size = 3.6, color = "black") +
    coord_flip(clip = "off") +
    scale_y_continuous(
      limits = c(0, upper),
      breaks = if (clamp01) seq(0, 1, 0.25) else waiver(),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey80"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, margin = margin(t = 2)),
      legend.position = "none",
      plot.margin = margin(10, 22, 10, 10)
    )

  ggsave(file.path(compare_root, "01_between_methods/plots", fname),
         p, width = 7.8, height = 4.9, dpi = 300)
}

save_bar(bm_df, "pearson",
         "Between-method Pearson correlation (replicate means: SSDM vs H2O)",
         "pearson_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "spearman",
         "Between-method Spearman correlation (replicate means: SSDM vs H2O)",
         "spearman_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "jaccard_hotspot_q",
         sprintf("Hotspot overlap (Jaccard, top %.0f%%, replicate means)", opts$q*100),
         "jaccard_q75_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "rmse",
         "Pixel-wise RMSE between methods (replicate means)",
         "rmse_bar.png", clamp01 = FALSE, digits = 3)

save_bar(bm_df, "mae",
         "Pixel-wise MAE between methods (replicate means)",
         "mae_bar.png", clamp01 = FALSE, digits = 3)

log_line(sprintf("Done. Outputs in: %s", compare_root), logf)
