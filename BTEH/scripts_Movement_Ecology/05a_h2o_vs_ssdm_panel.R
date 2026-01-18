# =============================================================================
# 05a_h2o_vs_ssdm_panel.R — Generate Before/After comparison panels
#
# Description: This script generates panels comparing habitat suitability maps
# before and after fence removal, including difference maps and GSL strips.
# =============================================================================


suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(gridExtra)
  library(grid)
  library(viridis)
})


# --- Robust sourcing ---------------------------------------------------------
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


source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))


`%||%` <- function(x, y) if (is.null(x)) y else x


# --- CLI ---------------------------------------------------------------------
opt <- list(
  make_option(c("--species"), type = "character", default = "E3,E4,E5", help = "Comma-separated: E3,E4,E5"),
  make_option(c("--methods"), type = "character", default = "H2O,SSDM", help = "Comma-separated: H2O,SSDM"),
  make_option(c("--thresholds"), type = "character", default = "0.25,0.5,0.75", help = "GSL thresholds"),
  make_option(c("--aoi"), type = "character", default = "data/shp/HV20233.shp", help = "Optional AOI vector; set '' to disable"),
  make_option(c("--mode"), type = "character", default = "FAST", help = "REPRO or FAST"),
  make_option(c("--metrics_csv"),
    type = "character", default = "results/compare/h2o_vs_ssdm/temporal_jaccard.csv",
    help = "Path to temporal_jaccard.csv"
  )
)
opts <- parse_args(OptionParser(option_list = opt))


cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)


species_list <- strsplit(opts$species, ",")[[1]]
methods_list <- strsplit(opts$methods, ",")[[1]]
thresholds_vec <- as.numeric(strsplit(opts$thresholds, ",")[[1]])


h2o_root <- cfg$paths$results_h2o %||% file.path("results", "H2O")
ssdm_root <- cfg$paths$results_ssdm %||% file.path("results", "SSDM")
compare_root <- file.path("results", "compare")
panels_root <- file.path("results", "panels")
logs_dir <- cfg$paths$logs %||% "logs"


dir_ensure(panels_root)
dir_ensure(logs_dir)
logf <- file.path(logs_dir, "06_before_after_panels.log")
log_line("Starting 06_before_after_panel_plots.R (using temporal_jaccard.csv)", logf)


# --- Load pre-computed metrics -----------------------------------------------
if (!file.exists(opts$metrics_csv)) {
  stop(sprintf("Metrics CSV not found: %s", opts$metrics_csv))
}
metrics_df <- read_csv(opts$metrics_csv, col_types = cols())
log_line(sprintf("Loaded metrics from: %s", opts$metrics_csv), logf)


# --- Load AOI if provided ----------------------------------------------------
aoi_vec <- NULL
if (!is.null(opts$aoi) && nchar(opts$aoi) > 0 && file.exists(opts$aoi)) {
  aoi_vec <- vect(opts$aoi)
  log_line(sprintf("AOI loaded: %s", opts$aoi), logf)
} else {
  log_line("No AOI crop (--aoi='' or file not found)", logf)
}


# --- Helper functions --------------------------------------------------------
align_rasters <- function(r1, r2) {
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- project(r2, crs(r1), method = "bilinear")
    r2 <- resample(r2, r1, method = "bilinear")
    r2 <- crop(r2, r1)
  }
  m <- !is.na(r1) & !is.na(r2)
  list(r1 = mask(r1, m, maskvalues = 0), r2 = mask(r2, m, maskvalues = 0))
}


crop_to_aoi <- function(r, aoi) {
  if (is.null(aoi)) {
    return(r)
  }
  aoi_proj <- project(aoi, crs(r))
  r <- crop(r, aoi_proj)
  r <- mask(r, aoi_proj)
  return(r)
}


plot_continuous_raster <- function(r, title, legend_title = "Suitability", use_viridis = TRUE) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "val")

  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = legend_title) +
    theme_void(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 8)),
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.width = unit(2.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )

  if (use_viridis) {
    p <- p + scale_fill_viridis_c(
      option = "viridis",
      na.value = NA,
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    )
  } else {
    # Centered diverging for difference maps
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p <- p + scale_fill_gradientn(
      colors = c("#0c97f3", "#ffffff", "#e41b1e"),
      limits = c(-lim, lim),
      na.value = "#ffffff",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    )
  }
  return(p)
}


plot_gsl_strip <- function(r_before, r_after, thr, title) {
  q_before <- quantile(values(r_before, mat = FALSE), thr, na.rm = TRUE)
  q_after <- quantile(values(r_after, mat = FALSE), thr, na.rm = TRUE)

  b_before <- r_before > q_before
  b_after <- r_after > q_after

  # 1=only_before(blue), 2=both(white), 3=only_after(red), 0=neither(gray)
  gsl <- ifel(
    b_before & b_after, 2,
    ifel(
      b_before & !b_after, 1,
      ifel(!b_before & b_after, 3, 0)
    )
  )

  df <- as.data.frame(gsl, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", "class")
  df <- df %>% filter(!is.na(class))
  df$class <- factor(df$class,
    levels = c(0, 1, 2, 3),
    labels = c("Neither", "Before only", "Both", "After only")
  )

  p <- ggplot(df, aes(x = x, y = y, fill = class)) +
    geom_raster() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    scale_fill_manual(
      values = c(
        "Neither" = "#f0f0f0", "Before only" = "#4575b4",
        "Both" = "#ffffff", "After only" = "#d73027"
      ),
      drop = FALSE
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 6)),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
  return(p)
}


make_metrics_table <- function(jaccard_vals, gain_vals, pct_delta_vals, thresholds_vec) {
  df <- data.frame(
    Metric = c(
      sprintf("J_t (thr=%.2f)", thresholds_vec[1]),
      sprintf("J_t (thr=%.2f)", thresholds_vec[2]),
      sprintf("J_t (thr=%.2f)", thresholds_vec[3]),
      sprintf("Gain km² (thr=%.2f)", thresholds_vec[1]),
      sprintf("Gain km² (thr=%.2f)", thresholds_vec[2]),
      sprintf("Gain km² (thr=%.2f)", thresholds_vec[3])
    ),
    Value = c(
      sprintf("%.3f", jaccard_vals[1]),
      sprintf("%.3f", jaccard_vals[2]),
      sprintf("%.3f", jaccard_vals[3]),
      sprintf("%.2f", gain_vals[1]),
      sprintf("%.2f", gain_vals[2]),
      sprintf("%.2f", gain_vals[3])
    )
  )

  # Color code rows: Jaccard 0-1 scale, Gain based on magnitude
  row_colors <- c(
    ifelse(jaccard_vals[1] > 0.5, "#a6d96a", ifelse(jaccard_vals[1] > 0.3, "#ffffbf", "#fdae61")),
    ifelse(jaccard_vals[2] > 0.5, "#a6d96a", ifelse(jaccard_vals[2] > 0.3, "#ffffbf", "#fdae61")),
    ifelse(jaccard_vals[3] > 0.5, "#a6d96a", ifelse(jaccard_vals[3] > 0.3, "#ffffbf", "#fdae61")),
    ifelse(gain_vals[1] > 30, "#a6d96a", ifelse(gain_vals[1] > 10, "#ffffbf", "#fdae61")),
    ifelse(gain_vals[2] > 30, "#a6d96a", ifelse(gain_vals[2] > 10, "#ffffbf", "#fdae61")),
    ifelse(gain_vals[3] > 30, "#a6d96a", ifelse(gain_vals[3] > 10, "#ffffbf", "#fdae61"))
  )

  tt <- ttheme_minimal(
    core = list(
      bg_params = list(fill = row_colors, col = "white", lwd = 2),
      fg_params = list(fontsize = 14, fontface = "bold", hjust = 0, x = 0.1)
    ),
    colhead = list(
      bg_params = list(fill = "gray30", col = "white", lwd = 2),
      fg_params = list(fontsize = 15, fontface = "bold", col = "white")
    ),
    rowhead = list(fg_params = list(fontsize = 0))
  )

  g <- tableGrob(df, rows = NULL, theme = tt)
  return(g)
}


# --- Main loop: species × methods --------------------------------------------
for (sp_base in species_list) {
  for (method in methods_list) {
    # CORRECTED: A = After, B = Before
    sp_after <- paste0(sp_base, "A") # After (fence removed)
    sp_before <- paste0(sp_base, "B") # Before (with fence)

    if (method == "H2O") {
      f_after <- file.path(h2o_root, "A", sp_after, sprintf("prediction_%s.tif", sp_after))
      f_before <- file.path(h2o_root, "B", sp_before, sprintf("prediction_%s.tif", sp_before))
      method_csv <- "h2o"
    } else if (method == "SSDM") {
      f_after <- file.path(ssdm_root, "A", sp_after, sprintf("ESDM_%s.tif", sp_after))
      f_before <- file.path(ssdm_root, "B", sp_before, sprintf("ESDM_%s.tif", sp_before))
      method_csv <- "SSDM"
    } else {
      next
    }

    if (!file.exists(f_after) || !file.exists(f_before)) {
      log_line(sprintf("Skip %s | %s — missing files", sp_base, method), logf)
      next
    }

    log_line(sprintf(
      "Processing: %s | %s (AOI crop: %s)",
      sp_base, method, ifelse(is.null(aoi_vec), "disabled", "enabled")
    ), logf)

    r_after <- rast(f_after)
    r_before <- rast(f_before)

    # AOI crop BEFORE alignment and all calculations
    r_after <- crop_to_aoi(r_after, aoi_vec)
    r_before <- crop_to_aoi(r_before, aoi_vec)

    al <- align_rasters(r_after, r_before)
    r_after <- al$r1
    r_before <- al$r2

    r_diff <- r_after - r_before

    # --- Row 1: Base maps + difference ---------------------------------------
    p1 <- plot_continuous_raster(r_before, sprintf("%s — %s (Before)", sp_base, method), "Suitability", TRUE)
    p2 <- plot_continuous_raster(r_after, sprintf("%s — %s (After)", sp_base, method), "Suitability", TRUE)
    p3 <- plot_continuous_raster(r_diff, "Δ (A − B)", "Δ Suitability", FALSE)

    # --- EXTRACT METRICS FROM CSV AND ALIGN WITH THRESHOLDS -------------------
    # Filter: ID == sp_base AND Method == method_csv
    sp_metrics <- metrics_df %>%
      filter(ID == sp_base, Method == method_csv)

    if (nrow(sp_metrics) != 3) {
      log_line(sprintf(
        "WARNING: Expected 3 rows for %s | %s, got %d. Skipping.",
        sp_base, method, nrow(sp_metrics)
      ), logf)
      next
    }

    # CRITICAL FIX: Match metrics to thresholds by row order in CSV
    # Assume CSV rows are in the SAME order as thresholds_vec (0.25, 0.5, 0.75)
    # If your CSV is not ordered by threshold, you need to explicitly sort/match

    jaccard_vals <- sp_metrics$J_t
    gain_vals <- sp_metrics$Gain_km2
    pct_delta_vals <- sp_metrics$pct_delta

    # Verify alignment: Print to log for debugging
    log_line(sprintf(
      "  Thresholds: [%.2f, %.2f, %.2f]",
      thresholds_vec[1], thresholds_vec[2], thresholds_vec[3]
    ), logf)
    log_line(sprintf(
      "  Temporal Jaccard: [%.3f, %.3f, %.3f]",
      jaccard_vals[1], jaccard_vals[2], jaccard_vals[3]
    ), logf)
    log_line(sprintf(
      "  Gain km²: [%.2f, %.2f, %.2f]",
      gain_vals[1], gain_vals[2], gain_vals[3]
    ), logf)

    # --- Row 2: GSL strips + metrics table -----------------------------------
    gsl_plots <- lapply(seq_along(thresholds_vec), function(i) {
      plot_gsl_strip(r_before, r_after, thresholds_vec[i], sprintf("GSL (thr=%.2f)", thresholds_vec[i]))
    })

    # --- Assemble layout (Table removed as per reviewer request) -------------
    row1 <- p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(1, 1, 1))
    row2 <- wrap_plots(gsl_plots[[1]], gsl_plots[[2]], gsl_plots[[3]],
      ncol = 3, widths = c(1, 1, 1)
    )

    final_panel <- row1 / row2 +
      plot_layout(heights = c(2, 1)) +
      plot_annotation(
        title = sprintf(
          "Elephant %s — %s | Before vs After: base maps, Δ (A−B), GSL strip, temporal metrics",
          sp_base, method
        ),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)))
      )

    # Save
    out_png <- file.path(panels_root, sprintf("%s_%s_before_after_panel.png", sp_base, method))
    ggsave(out_png, final_panel, width = 20, height = 12, dpi = 300, bg = "white")
    log_line(sprintf("Saved: %s", out_png), logf)
  }
}


log_line(sprintf("Done. Panels saved in: %s", panels_root), logf)
