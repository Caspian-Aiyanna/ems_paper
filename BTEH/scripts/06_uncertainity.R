#!/usr/bin/env Rscript
###############################################################################
# 06_uncertainity.R
# Paper-Ready Plotting (NO recomputation)
# Reads existing outputs under paper_results/ and produces grouped PNG figures.
#
# Assumptions (from your pipeline conventions):
# - Uncertainty outputs live under:
#     paper_results/uncertainty/<GROUP>/(H2O|SSDM)/
#     paper_results/uncertainty/<GROUP>/algo/
# - Before/After outputs live under:
#     paper_results/before_after/<ELEPHANT>/(H2O|SSDM)/
# - Pixelwise summaries & cross-group CSVs live under paper_results/uncertainty/ and
#   paper_results/before_after/
# - We only VISUALIZE existing .tif/.csv objects; we do not recompute anything.
###############################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
  library(tools)
})

# ---------------------- Config ----------------------
out_root    <- "paper_results"                             # where your existing outputs live
fig_root    <- file.path(out_root, "figures", "paper_ready")
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

# Groups/elephants used by your pipeline (edit if needed)
groups      <- c("E3A","E3B","E4A","E4B")
elephants   <- c("E3","E4")
methods     <- c("H2O","SSDM")
thr_main    <- 0.5

# ---------------------- Coastline (cached, reprojected on demand) ----------------------
.coast_cache <- new.env(parent = emptyenv())
get_coastline <- function(template_rast){
  key <- paste0("crs:", as.character(crs(template_rast)))
  if (!exists(key, envir = .coast_cache)){
    coast_sf <- tryCatch(
      rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf"),
      error = function(e) NULL
    )
    if (is.null(coast_sf)) return(NULL)
    coast_v  <- terra::vect(coast_sf)
    if (!is.na(crs(template_rast))) {
      coast_v <- tryCatch(terra::project(coast_v, template_rast), error = function(e) coast_v)
    }
    assign(key, coast_v, envir = .coast_cache)
  }
  get(key, envir = .coast_cache)
}

# ---------------------- Palette-aware map wrapper (no recompute) ----------------------
save_raster_map <- function(r, title, fname,
                            palette = c("sequential","diverging","binary"),
                            add_coast = TRUE, na_col = "grey85",
                            width_px = 1600, height_px = 1200, dpi = 240){
  palette <- match.arg(palette)
  dir.create(dirname(fname), TRUE, TRUE)

  v <- values(r); v <- v[is.finite(v)]
  png(filename = fname, width = width_px, height = height_px, res = dpi, bg = "white")
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mar = c(3.2, 3.6, 3.2, 6))

  coast <- if (add_coast) get_coastline(r) else NULL

  if (length(v) == 0){
    plot(r, main = title); if (!is.null(coast)) try(lines(crop(coast, ext(r)), lwd=0.4), silent=TRUE)
    return(invisible(NULL))
  }

  if (palette == "binary"){
    brks <- c(-0.5, 0.5, 1.5)
    cols <- c("#f0f0f0", "#3b528b")
    plot(r, main = title, col = cols, breaks = brks, axes = FALSE, box = TRUE, legend = TRUE, colNA = na_col)
  } else if (palette == "diverging"){
    rng  <- range(v, na.rm = TRUE); vmax <- max(abs(rng))
    cols <- grDevices::colorRampPalette(c("#3b528b", "#f7f7f7", "#b40426"))(201)
    plot(r, main = title, col = cols, zlim = c(-vmax, vmax), axes = FALSE, box = TRUE, legend = TRUE, colNA = na_col)
  } else {
    rng  <- range(v, na.rm = TRUE)
    cols <- viridisLite::viridis(200)
    plot(r, main = title, col = cols, zlim = rng, axes = FALSE, box = TRUE, legend = TRUE, colNA = na_col)
  }

  if (!is.null(coast)){
    try({ lines(crop(coast, ext(r)), lwd = 0.4, col = "black") }, silent = TRUE)
  }
}

# ---------------------- Utilities to *read* existing objects ----------------------
# Uncertainty (data-driven per method)
path_unc  <- function(g, m, f) file.path(out_root, "uncertainty", g, m, f)
# Uncertainty (algorithmic H2O vs SSDM)
path_algo <- function(g, f) file.path(out_root, "uncertainty", g, "algo", f)
# Before/After per elephant & method
path_ba   <- function(e, m, f) file.path(out_root, "before_after", e, m, f)

read_if_exists <- function(f){
  if (file.exists(f)) return(rast(f))
  stop("Missing file: ", f)
}

# ---------------------- Group panels: H2O vs SSDM maps per group ----------------------
plot_group_uncertainty_panel <- function(g){
  message("Group panel ", g)
  out_dir <- file.path(fig_root, "groups"); dir.create(out_dir, TRUE, TRUE)

  # Read rasters (stop loudly if any missing)
  r_mean_H <- read_if_exists(path_unc(g, "H2O",  "mean.tif"))
  r_sd_H   <- read_if_exists(path_unc(g, "H2O",  "sd.tif"))
  r_cv_H   <- read_if_exists(path_unc(g, "H2O",  "cv.tif"))
  r_st_H   <- read_if_exists(path_unc(g, "H2O",  sprintf("stability_thr%.2f.tif", thr_main)))

  r_mean_S <- read_if_exists(path_unc(g, "SSDM", "mean.tif"))
  r_sd_S   <- read_if_exists(path_unc(g, "SSDM", "sd.tif"))
  r_cv_S   <- read_if_exists(path_unc(g, "SSDM", "cv.tif"))
  r_st_S   <- read_if_exists(path_unc(g, "SSDM", sprintf("stability_thr%.2f.tif", thr_main)))

  # Panel 2 (cols) x 4 (rows): H2O vs SSDM for mean/sd/cv/stability
  fname <- file.path(out_dir, paste0(g, "_uncertainty_panel.png"))
  png(fname, width = 2200, height = 3600, res = 240, bg = "white")
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(4, 2), mar = c(3.2, 3.6, 3.2, 6))

  # Row 1: mean (sequential)
  plot(r_mean_H, main = paste0(g, " — H2O mean"),  col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_mean_H), ext(r_mean_H)), lwd=0.4), silent=TRUE)

  plot(r_mean_S, main = paste0(g, " — SSDM mean"), col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_mean_S), ext(r_mean_S)), lwd=0.4), silent=TRUE)

  # Row 2: sd (sequential)
  plot(r_sd_H,   main = paste0(g, " — H2O SD"),   col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_sd_H), ext(r_sd_H)), lwd=0.4), silent=TRUE)

  plot(r_sd_S,   main = paste0(g, " — SSDM SD"),  col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_sd_S), ext(r_sd_S)), lwd=0.4), silent=TRUE)

  # Row 3: cv (sequential)
  plot(r_cv_H,   main = paste0(g, " — H2O CV"),   col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_cv_H), ext(r_cv_H)), lwd=0.4), silent=TRUE)

  plot(r_cv_S,   main = paste0(g, " — SSDM CV"),  col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_cv_S), ext(r_cv_S)), lwd=0.4), silent=TRUE)

  # Row 4: stability (binary palette)
  plot(r_st_H,   main = paste0(g, " — H2O Stability (≥ ", thr_main, ")"),
       col=c("#f0f0f0", "#3b528b"), breaks=c(-0.5,0.5,1.5), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_st_H), ext(r_st_H)), lwd=0.4), silent=TRUE)

  plot(r_st_S,   main = paste0(g, " — SSDM Stability (≥ ", thr_main, ")"),
       col=c("#f0f0f0", "#3b528b"), breaks=c(-0.5,0.5,1.5), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_st_S), ext(r_st_S)), lwd=0.4), silent=TRUE)

  dev.off()
  message("Wrote: ", fname)
}

plot_group_algo_panel <- function(g){
  message("Algorithmic panel ", g)
  out_dir <- file.path(fig_root, "groups"); dir.create(out_dir, TRUE, TRUE)

  r_diff <- read_if_exists(path_algo(g, "diff_BminusA.tif"))         # SSDM − H2O (mean)
  r_abs  <- read_if_exists(path_algo(g, "absdiff.tif"))
  r_ui   <- read_if_exists(path_algo(g, "uncertainty_index_0to1.tif"))
  r_ag   <- read_if_exists(path_algo(g, sprintf("agreement_thr%.2f.tif", thr_main)))

  fname <- file.path(out_dir, paste0(g, "_algorithmic_panel.png"))
  png(fname, width = 2200, height = 2400, res = 240, bg = "white")
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(2, 2), mar = c(3.2,3.6,3.2,6))

  # Diff (diverging)
  v <- values(r_diff); v <- v[is.finite(v)]; vmax <- max(abs(range(v)))
  plot(r_diff, main = "SSDM − H2O (Mean)", col = colorRampPalette(c("#3b528b","#f7f7f7","#b40426"))(201),
       zlim=c(-vmax, vmax), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_diff), ext(r_diff)), lwd=0.4), silent=TRUE)

  # |Diff| (sequential)
  plot(r_abs,  main = "|Difference|", col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_abs), ext(r_abs)), lwd=0.4), silent=TRUE)

  # UI (sequential 0–1)
  plot(r_ui,   main = "Uncertainty Index (0–1)", col = viridisLite::viridis(200), zlim=c(0,1), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_ui), ext(r_ui)), lwd=0.4), silent=TRUE)

  # Agreement (binary)
  plot(r_ag,   main = paste0("Agreement (≥ ", thr_main, ")"),
       col=c("#f0f0f0", "#3b528b"), breaks=c(-0.5,0.5,1.5), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_ag), ext(r_ag)), lwd=0.4), silent=TRUE)

  dev.off(); message("Wrote: ", fname)
}

# ---------------------- Before/After 2×2 per elephant & method ----------------------
plot_before_after_panel <- function(e, m){
  message("Before/After panel ", e, " / ", m)
  out_dir <- file.path(fig_root, "before_after"); dir.create(out_dir, TRUE, TRUE)

  r_A   <- read_if_exists(path_ba(e, m, "A_mean.tif"))
  r_B   <- read_if_exists(path_ba(e, m, "B_mean.tif"))
  r_dif <- read_if_exists(path_ba(e, m, "BminusA_mean.tif"))

  # Change 3-class may not always be saved for non-median threshold; try 0.5
  f_change <- path_ba(e, m, sprintf("change_3class_thr%.2f.tif", thr_main))
  r_chg <- if (file.exists(f_change)) rast(f_change) else NULL

  fname <- file.path(out_dir, paste0(e, "_", m, "_A-B_panel.png"))
  png(fname, width = 2200, height = 2400, res = 240, bg = "white")
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(2, 2), mar = c(3.2,3.6,3.2,6))

  plot(r_A,   main = paste0(e, " / ", m, " — Before (A) mean"), col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_A), ext(r_A)), lwd=0.4), silent=TRUE)

  plot(r_B,   main = paste0(e, " / ", m, " — After (B) mean"),  col = viridisLite::viridis(200), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_B), ext(r_B)), lwd=0.4), silent=TRUE)

  # Diverging for difference
  v <- values(r_dif); v <- v[is.finite(v)]; vmax <- max(abs(range(v)))
  plot(r_dif, main = "After − Before (mean)", col = colorRampPalette(c("#3b528b","#f7f7f7","#b40426"))(201),
       zlim=c(-vmax, vmax), axes=FALSE, box=TRUE, legend=TRUE)
  try(lines(crop(get_coastline(r_dif), ext(r_dif)), lwd=0.4), silent=TRUE)

  if (!is.null(r_chg)){
    # 3 classes: loss(1)=#b40426, stable(2)=#f7f7f7, gain(3)=#3b528b
    plot(r_chg, main = paste0("Change classes (thr=", thr_main, ")"),
         col=c("#b40426", "#f7f7f7", "#3b528b"), axes=FALSE, box=TRUE, legend=TRUE)
    try(lines(crop(get_coastline(r_chg), ext(r_chg)), lwd=0.4), silent=TRUE)
  } else {
    plot(rast(r_A), main = "Change classes (missing)")
    mtext("change_3class raster not found", side=3, line=0.5, col="grey40")
  }

  dev.off(); message("Wrote: ", fname)
}

# ---------------------- Cross-group charts from existing CSVs ----------------------
plot_cross_group_charts <- function(){
  dir.create(file.path(fig_root, "cross_group"), TRUE, TRUE)
  out_dir <- file.path(fig_root, "cross_group")

  # Algorithmic global stats
  f_alg <- file.path(out_root, "uncertainty", "algorithmic_global_stats.csv")
  if (file.exists(f_alg)){
    alg <- read_csv(f_alg, show_col_types = FALSE)
    mk <- function(metric, ylab, fname){
      p <- ggplot(alg, aes(x=group, y=.data[[metric]]))+
        geom_col(fill="grey35")+
        labs(x="Group", y=ylab, title=paste(ylab, "by Group"))+
        theme_minimal(base_size = 13)
      ggsave(filename = file.path(out_dir, fname), plot = p, width = 8, height = 4.6, dpi = 300)
    }
    mk("pearson",       "Pearson (H2O vs SSDM)",     "alg_pearson_by_group.png")
    mk("rmse",          "RMSE (H2O vs SSDM)",        "alg_rmse_by_group.png")
    mk("mean_absdiff",  "Mean |difference|",         "alg_mean_absdiff_by_group.png")
    mk("kappa",         "Cohen's κ (thr=0.5)",       "alg_kappa_by_group.png")
    mk("iou",           "IoU (thr=0.5)",             "alg_iou_by_group.png")
  } else {
    message("algorithmic_global_stats.csv not found; skipping")
  }

  # Pixelwise summary
  f_pix <- file.path(out_root, "uncertainty", "pixelwise_summary_by_group_method.csv")
  if (file.exists(f_pix)){
    pix <- read_csv(f_pix, show_col_types = FALSE)
    p_cv <- ggplot(pix, aes(x=group, y=mean_cv, fill=method))+
      geom_col(position=position_dodge(0.7))+
      labs(x="Group", y="Mean CV (pixelwise)", title="Replicate Uncertainty (CV) by Group & Method")+
      theme_minimal(base_size = 13)
    ggsave(filename = file.path(out_dir, "replicate_cv_by_group_method.png"), plot = p_cv, width = 8.5, height = 4.6, dpi = 300)
  }

  # Before/After summaries
  f_jac  <- file.path(out_root, "before_after", "jaccard_all.csv")
  f_area <- file.path(out_root, "before_after", "area_change_all.csv")
  if (file.exists(f_jac)){
    jac <- read_csv(f_jac, show_col_types = FALSE)
    p_jac <- ggplot(jac, aes(x=factor(threshold), y=jaccard, fill=method))+
      geom_col(position=position_dodge(0.7))+
      facet_wrap(~ elephant, nrow=1)+
      labs(x="Threshold", y="Jaccard (A vs B)", title="Core Overlap (Before vs After)")+
      theme_minimal(base_size = 13)
    ggsave(filename = file.path(out_dir, "bva_jaccard_by_threshold.png"), plot = p_jac, width = 10, height = 4.8, dpi = 300)
  }
  if (file.exists(f_area)){
    area <- read_csv(f_area, show_col_types = FALSE)
    p_area <- ggplot(area, aes(x=factor(threshold), y=pct_change, fill=method))+
      geom_col(position=position_dodge(0.7))+
      facet_wrap(~ elephant, nrow=1)+
      labs(x="Threshold", y="% change (After vs Before)", title="Suitable Area Change")+
      theme_minimal(base_size = 13)
    ggsave(filename = file.path(out_dir, "bva_area_change_by_threshold.png"), plot = p_area, width = 10, height = 4.8, dpi = 300)
  }
}

# ---------------------- Drive it ----------------------
for (g in groups){
  # Per-group uncertainty + algorithmic panels
  try(plot_group_uncertainty_panel(g), silent = TRUE)
  try(plot_group_algo_panel(g), silent = TRUE)
}

for (e in elephants){
  for (m in methods){
    try(plot_before_after_panel(e, m), silent = TRUE)
  }
}

# Cross-group charts from CSVs you already wrote earlier
plot_cross_group_charts()

cat("\nAll paper-ready figure panels written under:\n", normalizePath(fig_root), "\n\n")


###############################################################################
# Paper-ready numbers & plots from existing outputs
# - Reads CSVs produced by your pipeline under paper_results/
# - Writes compact tables (CSV) and attractive figures (PNG)
#
# Outputs -> paper_results/figures/paper_ready_numbers/
#
# Requirements: ggplot2, dplyr, readr, tidyr, scales, viridisLite
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

# ---------------------- Paths ----------------------
root <- "paper_results"
out_dir <- file.path(root, "figures", "paper_ready_numbers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: clean theme
thm <- function(base=13){
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(color = "grey20"),
      legend.position = "right",
      legend.title = element_blank()
    )
}

# Palette helpers
pal_disc <- c("H2O" = "#2E86AB", "SSDM" = "#E07A5F")
pal_posneg <- c("Loss" = "#b40426", "Stable" = "#f7f7f7", "Gain" = "#3b528b")

# =============================================================================
# 1) AGREEMENT BETWEEN METHODS (H2O vs SSDM) — TABLE + HEATMAP
# =============================================================================
alg_csv <- file.path(root, "uncertainty", "algorithmic_global_stats.csv")
if (file.exists(alg_csv)) {
  alg <- read_csv(alg_csv, show_col_types = FALSE)
  # Compact paper-ready table
  tbl_alg <- alg %>%
    select(group, pearson, rmse, mean_absdiff, kappa, iou) %>%
    arrange(group)
  write_csv(tbl_alg, file.path(out_dir, "table_algorithmic_agreement.csv"))

  # Heatmap (values scaled column-wise 0–1 for visual comparability)
  alg_long <- tbl_alg %>%
    pivot_longer(-group, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>% mutate(value_scaled = rescale(value)) %>% ungroup()

  p_alg_heat <- ggplot(alg_long, aes(x = metric, y = group, fill = value_scaled))+
    geom_tile(color = "white", linewidth = 0.4)+
    geom_text(aes(label = ifelse(metric %in% c("pearson","kappa","iou"),
                                 scales::number(value, accuracy=0.01),
                                 scales::number(value, accuracy=0.001))),
              size = 6, color = "black")+
    scale_fill_gradientn(colours = viridis(200))+
    labs(title = "Algorithmic Agreement (H2O vs SSDM)",
         subtitle = "Numbers in tiles are raw values; fill = column-wise scaled",
         x = NULL, y = NULL, fill = "Scaled")+
    thm()
  ggsave(file.path(out_dir, "algorithmic_agreement_heatmap.png"), p_alg_heat, width = 10, height = 4.5, dpi = 300)
}

# =============================================================================
# 2) REPLICATE UNCERTAINTY (MEAN CV) — TABLE + GROUPED BAR
# =============================================================================
pix_csv <- file.path(root, "uncertainty", "pixelwise_summary_by_group_method.csv")
if (file.exists(pix_csv)) {
  pix <- read_csv(pix_csv, show_col_types = FALSE)
  # Parse group into elephant + period
  pix2 <- pix %>%
    mutate(elephant = substr(group, 1, 2), period = substr(group, 3, 3)) %>%
    mutate(period = ifelse(period == "A", "Before (A)", "After (B)"))

  # Paper-ready table (wide by period)
  tbl_cv <- pix2 %>%
    select(elephant, method, period, mean_cv) %>%
    mutate(mean_cv = round(mean_cv, 4)) %>%
    pivot_wider(names_from = period, values_from = mean_cv) %>%
    arrange(elephant, method)
  write_csv(tbl_cv, file.path(out_dir, "table_replicate_cv_by_elephant_method.csv"))

  # Grouped bars (A vs B by method within elephant)
  p_cv <- ggplot(pix2, aes(x = period, y = mean_cv, fill = method))+
    geom_col(position = position_dodge(width = 0.65), width = 0.6)+
    facet_wrap(~ elephant, nrow = 1)+
    scale_fill_manual(values = pal_disc)+
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05)))+
    labs(title = "Replicate Uncertainty (Mean CV)", y = "Mean CV (as %)", x = NULL)+
    thm()
  ggsave(file.path(out_dir, "replicate_cv_grouped_bars.png"), p_cv, width = 10, height = 4.2, dpi = 300)
}

# =============================================================================
# 3) BEFORE vs AFTER — AREA CHANGE + JACCARD + GAIN/LOSS/STABLE DECOMP
# =============================================================================
ba_dir <- file.path(root, "before_after")
jac_all_csv  <- file.path(ba_dir, "jaccard_all.csv")
area_all_csv <- file.path(ba_dir, "area_change_all.csv")

# ---- 3a. Line plot: Jaccard vs threshold ----
if (file.exists(jac_all_csv)) {
  jac <- read_csv(jac_all_csv, show_col_types = FALSE)
  p_jac <- ggplot(jac, aes(x = factor(threshold), y = jaccard, group = method, color = method))+
    geom_line(linewidth = 1.1, aes(group = method))+
    geom_point(size = 2.5)+
    facet_wrap(~ elephant, nrow = 1)+
    scale_color_manual(values = pal_disc)+
    scale_y_continuous(labels = percent_format(accuracy = 1))+
    labs(title = "Core Overlap (A vs B) across thresholds",
         x = "Threshold", y = "Jaccard (A ∩ B / A ∪ B)")+
    thm()
  ggsave(file.path(out_dir, "bva_jaccard_lineplot.png"), p_jac, width = 10, height = 4.2, dpi = 300)
}

# ---- 3b. Bars: % area change relative to A ----
if (file.exists(area_all_csv)) {
  area <- read_csv(area_all_csv, show_col_types = FALSE)
  p_area <- ggplot(area, aes(x = factor(threshold), y = pct_change, fill = method))+
    geom_col(position = position_dodge(width = 0.65), width = 0.6)+
    facet_wrap(~ elephant, nrow = 1)+
    scale_fill_manual(values = pal_disc)+
    scale_y_continuous(labels = function(x) paste0(x, "%"))+
    labs(title = "% Area Change (After vs Before)", x = "Threshold", y = "% change relative to A")+
    thm()
  ggsave(file.path(out_dir, "bva_area_change_bars.png"), p_area, width = 10, height = 4.2, dpi = 300)
}

# ---- 3c. Stacked bars: Loss / Stable / Gain as % of UNION ----
# We reconstruct intersection (I) and union (U) using A, B and J=I/U
# Given A, B and J: U = (A+B)/(1+J),  I = J*U
# Then: Loss = A - I ; Gain = B - I ; Stable = I
if (file.exists(area_all_csv) && file.exists(jac_all_csv)) {
  area <- read_csv(area_all_csv, show_col_types = FALSE)
  jac  <- read_csv(jac_all_csv,  show_col_types = FALSE)
  dat <- area %>% inner_join(jac, by = c("elephant","method","threshold")) %>%
    mutate(U = (area_A + area_B) / (1 + jaccard),
           I = jaccard * U,
           Loss = pmax(area_A - I, 0),
           Gain = pmax(area_B - I, 0),
           Stable = pmax(I, 0)) %>%
    mutate(total = Loss + Gain + Stable,
           Loss_pct = 100*Loss/total,
           Stable_pct = 100*Stable/total,
           Gain_pct = 100*Gain/total)

  # Save table at thr=0.5
  tbl_union <- dat %>% filter(abs(threshold - 0.5) < 1e-9) %>%
    transmute(elephant, method,
              Loss_pct = round(Loss_pct, 1), Stable_pct = round(Stable_pct, 1), Gain_pct = round(Gain_pct, 1)) %>%
    arrange(elephant, method)
  write_csv(tbl_union, file.path(out_dir, "table_gain_loss_stable_union_thr0.5.csv"))

  # Stacked bars across thresholds
  dat_long <- dat %>%
    select(elephant, method, threshold, Loss_pct, Stable_pct, Gain_pct) %>%
    pivot_longer(cols = ends_with("_pct"), names_to = "class", values_to = "pct") %>%
    mutate(class = recode(class, Loss_pct = "Loss", Stable_pct = "Stable", Gain_pct = "Gain"),
           class = factor(class, levels = c("Loss","Stable","Gain")))

  p_stack <- ggplot(dat_long, aes(x = factor(threshold), y = pct, fill = class))+
    geom_col(width = 0.65, color = "white", linewidth = 0.2)+
    facet_grid(elephant ~ method)+
    scale_fill_manual(values = pal_posneg)+
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,100))+
    labs(title = "Decomposition of Change (% of A ∪ B)", x = "Threshold", y = "% of union area")+
    thm()
  ggsave(file.path(out_dir, "bva_gain_stable_loss_stacked.png"), p_stack, width = 10.5, height = 6.6, dpi = 300)
}

# =============================================================================
# 4) COMPACT SUMMARY TABLE — ready to drop into manuscript
#    (one row per Elephant × Method at thr=0.5)
# =============================================================================
if (file.exists(area_all_csv) && file.exists(jac_all_csv) && file.exists(pix_csv)) {
  area <- read_csv(area_all_csv, show_col_types = FALSE)
  jac  <- read_csv(jac_all_csv,  show_col_types = FALSE)
  pix  <- read_csv(pix_csv,      show_col_types = FALSE)

  # mean CV by elephant/method and period (A,B)
  cv_tbl <- pix %>% mutate(elephant = substr(group,1,2), period = substr(group,3,3)) %>%
    select(elephant, method, period, mean_cv) %>%
    mutate(period = ifelse(period=="A","CV_A","CV_B")) %>%
    pivot_wider(names_from = period, values_from = mean_cv, values_fill = NA_real_)

  # derive union-based components at thr=0.5
  comp <- area %>% inner_join(jac, by = c("elephant","method","threshold")) %>%
    filter(abs(threshold - 0.5) < 1e-9) %>%
    mutate(U=(area_A+area_B)/(1+jaccard), I=jaccard*U,
           Loss=area_A-I, Gain=area_B-I, Stable=I,
           Loss_pct=100*Loss/(Loss+Gain+Stable),
           Stable_pct=100*Stable/(Loss+Gain+Stable),
           Gain_pct=100*Gain/(Loss+Gain+Stable)) %>%
    select(elephant, method, pct_change, Loss_pct, Stable_pct, Gain_pct)

  final_tbl <- comp %>% left_join(cv_tbl, by = c("elephant","method")) %>%
    transmute(Elephant = elephant, Method = method,
              `Δ area % (B vs A)` = round(pct_change, 1),
              `Loss % (union)` = round(Loss_pct, 1),
              `Stable % (union)` = round(Stable_pct, 1),
              `Gain % (union)` = round(Gain_pct, 1),
              `CV_A` = round(CV_A, 3), `CV_B` = round(CV_B, 3)) %>%
    arrange(Elephant, Method)

  write_csv(final_tbl, file.path(out_dir, "table_compact_summary_thr0.5.csv"))
}

cat("\nAll numeric tables and figures written to:\n", normalizePath(out_dir), "\n\n")
