# =============================================================================
# 05_h2o_vs_ssdm.R — Between-method comparison (H2O vs SSDM)
# - Uses config.yml + utils
# - Per-run A/B and optional single species
# - Outputs: CSV metrics, diff rasters, maps, annotated bar charts
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
  library(grid)    # unit() for legend guide
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

# --- CLI ---
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

# --- paths & logs ---
h2o_root   <- cfg$paths$results_h2o %||% file.path("results","H2O")
ssdm_root  <- cfg$paths$results_ssdm %||% file.path("results","SSDM")
compare_root <- file.path("results","compare", toupper(opts$run))
logs_dir   <- cfg$paths$logs %||% "logs"
logf       <- file.path(logs_dir, sprintf("05_compare_%s.log", toupper(opts$run)))
dir_ensure(compare_root); dir_ensure(logs_dir)

dir_ensure(file.path(compare_root, "01_between_methods/rasters"))
dir_ensure(file.path(compare_root, "01_between_methods/plots"))
dir_ensure(file.path(compare_root, "03_maps"))

log_line(sprintf("Starting 05_h2o_vs_ssdm.R  (mode=%s, run=%s, q=%.2f)",
                 cfg$mode, toupper(opts$run), opts$q), logf)

# --- discovery of species present in both methods ---
# H2O: results/H2O/<RUN>/<SP>/prediction_<SP>.tif
h2o_dir_run <- file.path(h2o_root, toupper(opts$run))
h2o_sp <- list.dirs(h2o_dir_run, full.names = FALSE, recursive = FALSE)
h2o_pred_paths <- setNames(
  file.path(h2o_dir_run, h2o_sp, paste0("prediction_", h2o_sp, ".tif")),
  h2o_sp
)

# SSDM: results/SSDM/Outputs/<SP>/ESDM_<SP>.tif
ssdm_out_dir <- file.path(ssdm_root, "Outputs")
ssdm_sp <- list.dirs(ssdm_out_dir, full.names = FALSE, recursive = FALSE)
ssdm_pred_paths <- setNames(
  file.path(ssdm_out_dir, ssdm_sp, paste0("ESDM_", ssdm_sp, ".tif")),
  ssdm_sp
)

# Species that have both rasters present
both_sp <- intersect(names(h2o_pred_paths)[file.exists(h2o_pred_paths)],
                     names(ssdm_pred_paths)[file.exists(ssdm_pred_paths)])
if (!is.null(opts$species)) both_sp <- intersect(both_sp, opts$species)
stopifnot(length(both_sp) > 0)
log_line(sprintf("Species compared: %s", paste(both_sp, collapse = ", ")), logf)

# --- helpers ------------------------------------------------------------------
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
  barheight_cm = 10,
  barwidth_cm  = 0.6,
  base_size    = 12,
  title_face   = "bold"
) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
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

# --- main loop: compute metrics, write maps + diffs ---------------------------
bm_rows <- list()

for (sp in both_sp) {
  f_h2o  <- h2o_pred_paths[[sp]]
  f_ssdm <- ssdm_pred_paths[[sp]]
  log_line(sprintf("Compare %s | H2O: %s | SSDM: %s", sp, basename(f_h2o), basename(f_ssdm)), logf)

  r_h2o  <- rast(f_h2o)
  r_ssdm <- rast(f_ssdm)
  al     <- align_to(r_ssdm, r_h2o, categorical = FALSE)
  r_ssdm <- al$r1; r_h2o <- al$r2

  # Save maps with larger legends
  plot_raster_continuous(r_ssdm, paste0(sp, " — SSDM suitability"),
                         file.path(compare_root, "03_maps", paste0(sp, "_SSDM.png")))
  plot_raster_continuous(r_h2o,  paste0(sp, " — H2O suitability"),
                         file.path(compare_root, "03_maps", paste0(sp, "_H2O.png")))

  # Metrics on valid pixels
  v1 <- values(r_ssdm, mat = FALSE)
  v2 <- values(r_h2o,  mat = FALSE)
  idx <- !is.na(v1) & !is.na(v2)
  pearson  <- suppressWarnings(cor(v1[idx], v2[idx], method = "pearson"))
  spearman <- suppressWarnings(cor(v1[idx], v2[idx], method = "spearman"))
  rmse     <- sqrt(mean((v1[idx] - v2[idx])^2))
  mae      <- mean(abs(v1[idx] - v2[idx]))

  # Hotspot Jaccard at quantile q
  q1 <- quantile(v1[idx], opts$q, na.rm = TRUE)
  q2 <- quantile(v2[idx], opts$q, na.rm = TRUE)
  h1 <- r_ssdm > q1
  h2 <- r_h2o  > q2
  jac_hot <- jaccard_binary(h1, h2)

  # Diff raster (H2O − SSDM)
  diff_r <- r_h2o - r_ssdm
  writeRaster(diff_r,
              file.path(compare_root, "01_between_methods/rasters",
                        paste0("diff_H2O_minus_SSDM_", sp, ".tif")),
              overwrite = TRUE,
              gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"),
              datatype = "FLT4S")
  plot_raster_continuous(diff_r,
                         paste0("Between-method difference (H2O − SSDM): ", sp),
                         file.path(compare_root, "01_between_methods/plots",
                                   paste0("diff_H2O_minus_SSDM_", sp, ".png")),
                         center0 = TRUE,
                         legend_title = "Δ suitability")

  bm_rows[[length(bm_rows) + 1L]] <- tibble(
    dataset = sp,
    pearson = pearson, spearman = spearman,
    rmse = rmse, mae = mae, jaccard_hotspot_q = jac_hot
  )
}

bm_df <- bind_rows(bm_rows)
readr::write_csv(bm_df,
                 file.path(compare_root, "01_between_methods", "between_method_metrics.csv"))

# --- annotated bar charts -----------------------------------------------------
mpl_blue <- "#1f77b4"

# Order that looks nice; only keep present species
target_order <- c("E1B","E5B","E5A","E2B","E4A","E3B","E3A","E4B","E6B")
present_order <- intersect(target_order, bm_df$dataset)
bm_df <- bm_df %>% mutate(dataset = factor(dataset, levels = rev(present_order))) %>% arrange(dataset)

save_bar <- function(df, metric, title, fname, clamp01 = TRUE, digits = 3) {
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

  ggsave(file.path(compare_root, "01_between_methods/plots", fname),
         p, width = 7.6, height = 4.8, dpi = 300)
}

save_bar(bm_df, "pearson",
         "Between-method Pearson correlation (SSDM vs H2O)",
         "pearson_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "spearman",
         "Between-method Spearman correlation (SSDM vs H2O)",
         "spearman_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "jaccard_hotspot_q",
         sprintf("Hotspot overlap (Jaccard, top %.0f%%)", opts$q*100),
         "jaccard_q_bar.png", clamp01 = TRUE, digits = 3)

save_bar(bm_df, "rmse",
         "Pixel-wise RMSE between methods",
         "rmse_bar.png", clamp01 = FALSE, digits = 3)

save_bar(bm_df, "mae",
         "Pixel-wise MAE between methods",
         "mae_bar.png", clamp01 = FALSE, digits = 3)

log_line(sprintf("Done. Outputs in: %s", compare_root), logf)
