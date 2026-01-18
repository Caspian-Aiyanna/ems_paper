#!/usr/bin/env Rscript
# =============================================================================
# 05_h2o_vs_ssdm.R — Between-method comparison (H2O vs SSDM)
# - Config + utils driven
# - If --run is omitted, auto-runs BOTH A and B by self-invocation
# - Optional --species to limit to one dataset (e.g., E3A)
# - Inputs:
#     H2O : results/H2O/<RUN>/<SP>/prediction_<SP>.tif
#     SSDM: results/SSDM/<RUN>/<SP>/ESDM_<SP>.tif
# - Outputs (all three locations are written):
#     Central: results/compare/h2o_vs_ssdm/<RUN>/
#     H2O per-SP:  results/H2O/<RUN>/<SP>/compare/
#     SSDM per-SP: results/SSDM/<RUN>/<SP>/compare/
#   (metrics CSV, diff rasters, hotspot rasters, maps, per-elephant panels)
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
  library(grid) # unit() for legend guide
  library(patchwork) # panels
  library(sf)
})

# --------------------------- robust script path -------------------------------
.this_file <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  filearg <- grep("^--file=", ca, value = TRUE)
  if (length(filearg)) {
    return(normalizePath(sub("^--file=", "", filearg), mustWork = FALSE))
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

# --------------------------- project utilities --------------------------------
source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# --------------------------- CLI ----------------------------------------------
opt_list <- list(
  make_option(c("--run"), type = "character", default = NULL, help = "Run tag: A or B. If omitted, the script runs BOTH A and B."),
  make_option(c("--mode"), type = "character", default = NULL, help = "REPRO or FAST (overrides config.yml)"),
  make_option(c("--species"), type = "character", default = NULL, help = "Optional single dataset (e.g., E3A)"),
  make_option(c("--q"), type = "double", default = 0.75, help = "Hotspot quantile (default 0.75)"),
  make_option(c("--aoi"), type = "character", default = "data/shp/HV20233.shp", help = "Optional AOI vector; set '' to disable")
)
opts <- parse_args(OptionParser(option_list = opt_list))

# --- Self-invocation when --run is omitted ------------------------------------
if (is.null(opts$run)) {
  message("[INFO] No --run provided; running both A and B sequentially.")

  script <- .this_file()
  if (!nzchar(script) || is.na(script)) stop("[FATAL] Cannot resolve script path for self-invocation.")

  # rebuild trailing args, dropping any existing --run/-r pair if present
  ta <- commandArgs(trailingOnly = TRUE)
  drop_idx <- which(ta %in% c("--run", "-r"))
  if (length(drop_idx)) {
    keep <- setdiff(seq_along(ta), union(drop_idx, drop_idx + 1L))
    ta <- ta[keep]
  }

  for (r in c("A", "B")) {
    message(sprintf("[INFO] -> Launching child for run %s ...", r))
    system2(file.path(R.home("bin"), "Rscript"),
      args = c(script, ta, "--run", r),
      stdout = "", stderr = ""
    )
  }
  quit(status = 0)
}

# --------------------------- Config & mode ------------------------------------
opts$run <- toupper(opts$run)
if (!is.null(opts$mode)) opts$mode <- toupper(opts$mode)

cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- opts$mode
mode <- set_mode(cfg) # sets seeds/cores per REPRO|FAST

# --------------------------- Paths & logging ----------------------------------
h2o_root <- cfg$paths$results_h2o %||% file.path(.root, "results", "H2O")
ssdm_root <- cfg$paths$results_ssdm %||% file.path(.root, "results", "SSDM")
compare_root <- file.path(.root, "results", "compare", "h2o_vs_ssdm", opts$run)
logs_dir <- cfg$paths$logs %||% file.path(.root, "logs")
dir_ensure(compare_root)
dir_ensure(logs_dir)

# central subdirs
dir_ensure(file.path(compare_root, "01_between_methods/rasters"))
dir_ensure(file.path(compare_root, "01_between_methods/plots"))
dir_ensure(file.path(compare_root, "02_hotspots"))
dir_ensure(file.path(compare_root, "03_maps"))
dir_ensure(file.path(compare_root, "04_panels"))

logf <- file.path(logs_dir, sprintf("05_compare_%s.log", opts$run))
log_line(sprintf("Starting 05_h2o_vs_ssdm.R  (mode=%s, run=%s, q=%.2f)", cfg$mode, opts$run, opts$q), logf)

# --------------------------- Discover datasets --------------------------------
# H2O: results/H2O/<RUN>/<SP>/prediction_<SP>.tif
h2o_dir_run <- file.path(h2o_root, opts$run)
if (!dir.exists(h2o_dir_run)) {
  stop(sprintf("[FATAL] H2O run directory does not exist: %s", h2o_dir_run))
}
h2o_sp <- list.dirs(h2o_dir_run, full.names = FALSE, recursive = FALSE)
h2o_pred_paths <- setNames(file.path(h2o_dir_run, h2o_sp, paste0("prediction_", h2o_sp, ".tif")), h2o_sp)

# SSDM: results/SSDM/<RUN>/<SP>/ESDM_<SP>.tif
ssdm_out_dir <- file.path(ssdm_root, opts$run)
if (!dir.exists(ssdm_out_dir)) {
  stop(sprintf("[FATAL] SSDM run directory does not exist: %s", ssdm_out_dir))
}
ssdm_sp <- list.dirs(ssdm_out_dir, full.names = FALSE, recursive = FALSE)
ssdm_pred_paths <- setNames(file.path(ssdm_out_dir, ssdm_sp, paste0("ESDM_", ssdm_sp, ".tif")), ssdm_sp)

both_sp <- intersect(
  names(h2o_pred_paths)[file.exists(h2o_pred_paths)],
  names(ssdm_pred_paths)[file.exists(ssdm_pred_paths)]
)

if (!is.null(opts$species)) {
  # allow single or comma-separated
  want <- toupper(unlist(strsplit(opts$species, "[,;\\s]+")))
  both_sp <- intersect(both_sp, want)
}
stopifnot(length(both_sp) > 0)
log_line(sprintf("Species compared: %s", paste(both_sp, collapse = ", ")), logf)

# --------------------------- Optional AOI -------------------------------------
aoi <- NULL
if (!is.null(opts$aoi) && nzchar(opts$aoi) && file.exists(opts$aoi)) {
  aoi <- tryCatch(suppressWarnings(sf::st_read(opts$aoi, quiet = TRUE)), error = function(e) NULL)
  if (!is.null(aoi)) log_line(sprintf("[AOI] Loaded %s", opts$aoi), logf)
}

# --------------------------- Helpers ------------------------------------------
align_to <- function(r1, r2, categorical = FALSE) {
  # align r2 to r1 (CRS, res, extent), then mask to common NA
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
  barwidth_cm = 0.7,
  base_size = 12,
  title_face = "bold",
  aoi = NULL
) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "val")
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = legend_title) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, face = title_face),
      legend.position = "right",
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size)
    )

  if (center0) {
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p <- p + scale_fill_gradientn(
      colors = c("#2c7fb8", "#ffffbf", "#d7191c"),
      limits = c(-lim, lim),
      na.value = NA,
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth = unit(barwidth_cm, "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  } else {
    p <- p + scale_fill_viridis_c(
      na.value = NA,
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth = unit(barwidth_cm, "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  }

  # AOI overlay (safe; no aes inheritance)
  if (!is.null(aoi)) {
    p <- p + geom_sf(data = aoi, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)
  }

  ggsave(out_png, p, width = 8.6, height = 6.4, dpi = 300)
}

save_bar <- function(df, metric, title, fname, outdir, clamp01 = TRUE, digits = 3) {
  mpl_blue <- "#1f77b4"
  df <- df %>% mutate(label = sprintf(paste0("%.", digits, "f"), .data[[metric]]))
  maxv <- suppressWarnings(max(df[[metric]], na.rm = TRUE))
  upper <- if (clamp01) 1 else max(1e-9, maxv) * 1.10

  p <- ggplot(df, aes(x = dataset, y = .data[[metric]])) +
    geom_col(fill = mpl_blue, width = 0.8) +
    geom_text(aes(label = label), hjust = -0.15, size = 3.4, color = "black") +
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
      plot.margin = margin(10, 20, 10, 10)
    )

  dir_ensure(outdir)
  ggsave(file.path(outdir, fname), p, width = 7.6, height = 4.8, dpi = 300)
}

# --------------------------- Main loop ----------------------------------------
bm_rows <- list()

for (sp in both_sp) {
  f_h2o <- h2o_pred_paths[[sp]]
  f_ssdm <- ssdm_pred_paths[[sp]]
  log_line(sprintf("Compare %s | H2O: %s | SSDM: %s", sp, basename(f_h2o), basename(f_ssdm)), logf)

  # per-species compare folders
  sp_h2o_cmp <- file.path(h2o_root, opts$run, sp, "compare")
  sp_ssdm_cmp <- file.path(ssdm_root, opts$run, sp, "compare")
  dir_ensure(sp_h2o_cmp)
  dir_ensure(sp_ssdm_cmp)

  r_h2o <- rast(f_h2o)
  r_ssdm <- rast(f_ssdm)
  al <- align_to(r_ssdm, r_h2o, categorical = FALSE)
  r_ssdm <- al$r1
  r_h2o <- al$r2

  # --- maps with larger legends
  map_ssdm <- file.path(compare_root, "03_maps", paste0(sp, "_SSDM.png"))
  map_h2o <- file.path(compare_root, "03_maps", paste0(sp, "_H2O.png"))
  plot_raster_continuous(r_ssdm, paste0(sp, " — SSDM suitability"), map_ssdm, aoi = aoi)
  plot_raster_continuous(r_h2o, paste0(sp, " — H2O suitability"), map_h2o, aoi = aoi)

  # copy maps to per-species folders
  file.copy(map_ssdm, file.path(sp_h2o_cmp, basename(map_ssdm)), overwrite = TRUE)
  file.copy(map_ssdm, file.path(sp_ssdm_cmp, basename(map_ssdm)), overwrite = TRUE)
  file.copy(map_h2o, file.path(sp_h2o_cmp, basename(map_h2o)), overwrite = TRUE)
  file.copy(map_h2o, file.path(sp_ssdm_cmp, basename(map_h2o)), overwrite = TRUE)

  # --- metrics on valid pixels
  v1 <- values(r_ssdm, mat = FALSE)
  v2 <- values(r_h2o, mat = FALSE)
  idx <- !is.na(v1) & !is.na(v2)
  pearson <- suppressWarnings(cor(v1[idx], v2[idx], method = "pearson"))
  spearman <- suppressWarnings(cor(v1[idx], v2[idx], method = "spearman"))
  rmse <- sqrt(mean((v1[idx] - v2[idx])^2))
  mae <- mean(abs(v1[idx] - v2[idx]))

  # --- Hotspot Jaccard at quantile q
  q1 <- quantile(v1[idx], opts$q, na.rm = TRUE)
  q2 <- quantile(v2[idx], opts$q, na.rm = TRUE)
  h1 <- r_ssdm > q1
  h2 <- r_h2o > q2
  jac_hot <- jaccard_binary(h1, h2)

  # write hotspot rasters (central + per-species)
  hs_ssdm <- file.path(compare_root, "02_hotspots", paste0("hotspot_ssdm_q", round(opts$q * 100), "_", sp, ".tif"))
  hs_h2o <- file.path(compare_root, "02_hotspots", paste0("hotspot_h2o_q", round(opts$q * 100), "_", sp, ".tif"))
  hs_ovl <- file.path(compare_root, "02_hotspots", paste0("hotspot_overlap_q", round(opts$q * 100), "_", sp, ".tif"))

  writeRaster(h1, hs_ssdm, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(h2, hs_h2o, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(h1 & h2, hs_ovl, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  file.copy(hs_ssdm, file.path(sp_h2o_cmp, basename(hs_ssdm)), overwrite = TRUE)
  file.copy(hs_ssdm, file.path(sp_ssdm_cmp, basename(hs_ssdm)), overwrite = TRUE)
  file.copy(hs_h2o, file.path(sp_h2o_cmp, basename(hs_h2o)), overwrite = TRUE)
  file.copy(hs_h2o, file.path(sp_ssdm_cmp, basename(hs_h2o)), overwrite = TRUE)
  file.copy(hs_ovl, file.path(sp_h2o_cmp, basename(hs_ovl)), overwrite = TRUE)
  file.copy(hs_ovl, file.path(sp_ssdm_cmp, basename(hs_ovl)), overwrite = TRUE)

  # --- Diff raster (H2O − SSDM)
  diff_r <- r_h2o - r_ssdm
  diff_tif <- file.path(compare_root, "01_between_methods/rasters", paste0("diff_H2O_minus_SSDM_", sp, ".tif"))
  writeRaster(diff_r, diff_tif,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"),
    datatype = "FLT4S"
  )

  diff_png <- file.path(compare_root, "01_between_methods/plots", paste0("diff_H2O_minus_SSDM_", sp, ".png"))
  plot_raster_continuous(diff_r,
    paste0("Between-method difference (H2O − SSDM): ", sp),
    diff_png,
    center0 = TRUE,
    legend_title = "Δ suitability",
    aoi = aoi
  )

  # copy diff to per-species folders
  file.copy(diff_tif, file.path(sp_h2o_cmp, basename(diff_tif)), overwrite = TRUE)
  file.copy(diff_tif, file.path(sp_ssdm_cmp, basename(diff_tif)), overwrite = TRUE)
  file.copy(diff_png, file.path(sp_h2o_cmp, basename(diff_png)), overwrite = TRUE)
  file.copy(diff_png, file.path(sp_ssdm_cmp, basename(diff_png)), overwrite = TRUE)

  # --- Accumulate row
  bm_rows[[length(bm_rows) + 1L]] <- tibble(
    Dataset = sp,
    `Pearson r` = pearson,
    `Spearman ρ` = spearman,
    RMSE = rmse,
    MAE = mae,
    `Jaccard_0.75` = jac_hot
  )

  # --- Per-elephant PANEL -----------------------------------------------------
  as_xy <- function(r) {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(df) <- c("x", "y", "val")
    df
  }
  df_ssdm <- as_xy(r_ssdm)
  df_h2o <- as_xy(r_h2o)
  df_diff <- as_xy(diff_r)
  df_ag <- as_xy((h1 & h2) * 1)

  p_ssdm <- ggplot(df_ssdm, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(title = "SSDM suitability", fill = "Suitability") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())
  if (!is.null(aoi)) p_ssdm <- p_ssdm + geom_sf(data = aoi, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)

  p_h2o <- ggplot(df_h2o, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(title = "H2O suitability", fill = "Suitability") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())
  if (!is.null(aoi)) p_h2o <- p_h2o + geom_sf(data = aoi, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)

  lim <- max(abs(range(df_diff$val, na.rm = TRUE)))
  p_diff <- ggplot(df_diff, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    scale_fill_gradientn(
      colors = c("#2c7fb8", "#ffffbf", "#d7191c"),
      limits = c(-lim, lim)
    ) +
    labs(title = "Δ (H2O − SSDM)", fill = "Δ") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())
  if (!is.null(aoi)) p_diff <- p_diff + geom_sf(data = aoi, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)

  pct_agree <- round(100 * mean(df_ag$val == 1, na.rm = TRUE), 2)
  p_agree <- ggplot(df_ag, aes(x = x, y = y, fill = factor(val))) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(values = c("white", "grey50"), breaks = c("0", "1"), labels = c("No", "Yes")) +
    labs(
      title = paste0("Hotspot Overlap @", round(opts$q * 100), "% (", pct_agree, "%)"),
      fill = "Overlap"
    ) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())
  if (!is.null(aoi)) p_agree <- p_agree + geom_sf(data = aoi, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)

  mtxt <- paste0(
    "Pearson: ", sprintf("%.3f", pearson), "\n",
    "Spearman: ", sprintf("%.3f", spearman), "\n",
    "RMSE: ", sprintf("%.3f", rmse), "\n",
    "MAE: ", sprintf("%.3f", mae), "\n",
    "Jaccard@", round(opts$q * 100), "%: ", sprintf("%.3f", jac_hot)
  )
  p_metrics <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, label = mtxt, size = 4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    labs(title = "Metrics")

  panel <- (p_ssdm | p_h2o) / (p_diff | p_agree) / p_metrics +
    plot_annotation(title = paste0("H2O vs SSDM — ", sp, " (Run ", opts$run, ")"))

  panel_out_central <- file.path(compare_root, "04_panels", paste0("panel_", sp, ".png"))
  ggsave(panel_out_central, panel, width = 12, height = 12, dpi = 300)

  # copy panels to per-species
  file.copy(panel_out_central, file.path(sp_h2o_cmp, basename(panel_out_central)), overwrite = TRUE)
  file.copy(panel_out_central, file.path(sp_ssdm_cmp, basename(panel_out_central)), overwrite = TRUE)
}

# --------------------------- Metrics CSV + bars -------------------------------
bm_df <- bind_rows(bm_rows) %>%
  arrange(Dataset) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

# Write CSV with proper column names matching LaTeX table
metrics_csv <- file.path(compare_root, "01_between_methods", "between_method_metrics.csv")
readr::write_csv(bm_df, metrics_csv)
log_line(sprintf("Metrics CSV written: %s", metrics_csv), logf)

# Copy to per-species folders
for (sp in bm_df$Dataset) {
  sp_h2o_cmp <- file.path(h2o_root, opts$run, sp, "compare")
  sp_ssdm_cmp <- file.path(ssdm_root, opts$run, sp, "compare")
  file.copy(metrics_csv, file.path(sp_h2o_cmp, "between_method_metrics.csv"), overwrite = TRUE)
  file.copy(metrics_csv, file.path(sp_ssdm_cmp, "between_method_metrics.csv"), overwrite = TRUE)
}

# --- Generate 5 separate metric plots for all datasets -----------------------
bar_dir <- file.path(compare_root, "01_between_methods/plots")

# Prepare data for plotting (need internal column names)
plot_df <- bm_df %>%
  rename(
    dataset = Dataset,
    pearson = `Pearson r`,
    spearman = `Spearman ρ`,
    rmse = RMSE,
    mae = MAE,
    jaccard_hotspot_q = `Jaccard_0.75`
  )

save_bar(plot_df, "pearson",
  "Between-method Pearson correlation (SSDM vs H2O)",
  "pearson_bar.png", bar_dir,
  clamp01 = TRUE, digits = 3
)

save_bar(plot_df, "spearman",
  "Between-method Spearman correlation (SSDM vs H2O)",
  "spearman_bar.png", bar_dir,
  clamp01 = TRUE, digits = 3
)

save_bar(plot_df, "jaccard_hotspot_q",
  sprintf("Hotspot overlap (Jaccard, top %.0f%%)", opts$q * 100),
  "jaccard_bar.png", bar_dir,
  clamp01 = TRUE, digits = 3
)

save_bar(plot_df, "rmse",
  "Pixel-wise RMSE between methods",
  "rmse_bar.png", bar_dir,
  clamp01 = FALSE, digits = 3
)

save_bar(plot_df, "mae",
  "Pixel-wise MAE between methods",
  "mae_bar.png", bar_dir,
  clamp01 = FALSE, digits = 3
)

log_line(sprintf("Done. Central outputs in: %s", compare_root), logf)
log_line("All 5 metric plots generated successfully", logf)
log_line("=========================================", logf)
