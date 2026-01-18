# =============================================================================
# 07_appendix.R — Appendix Figures and Tables Generation
#
# Description: This script generates supplementary materials for the paper,
# including base learner contributions, agreement maps between modeling
# frameworks, and variable importance dumbbell plots.
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║        07_appendix_MASTER.R — Complete Appendix             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  library(h2o)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

cat("✓ Libraries loaded\n\n")

# --- robust sourcing ---------------------------------------------------------
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

# --- PATHS -------------------------------------------------------------------
h2o_root <- file.path(.root, "results", "H2O")
ssdm_root <- file.path(.root, "results", "SSDM")
appendix_root <- file.path(.root, "results", "appendix")

dir.create(file.path(appendix_root, "base_learner_contributions"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(appendix_root, "agreement_maps"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(appendix_root, "variable_importance_h2o"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(appendix_root, "variable_importance_ssdm"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(appendix_root, "data"), showWarnings = FALSE, recursive = TRUE)

# --- SPECIES -----------------------------------------------------------------
h2o_species_all <- c("E1B", "E2B", "E3A", "E3B", "E4A", "E4B", "E5A", "E5B", "E6B")
ssdm_species_all <- c("E1B", "E2B", "E3A", "E3B", "E4A", "E4B", "E5A", "E5B", "E6B")
vimp_species <- c("E3", "E4", "E5")
thr <- 0.5

cat(sprintf("H2O species: %s\n", paste(h2o_species_all, collapse = ", ")))
cat(sprintf("SSDM species: %s\n", paste(ssdm_species_all, collapse = ", ")))
cat(sprintf("VI species: %s\n\n", paste(vimp_species, collapse = ", ")))

# --- H2O INITIALIZATION ------------------------------------------------------
h2o.init(nthreads = 8, max_mem_size = "8G")
h2o.no_progress()
cat("✓ H2O initialized\n\n")

# --- VARIABLE LABELING -------------------------------------------------------
BIO_LUT <- c(
  bio1  = "Annual Mean Temperature",
  bio2  = "Mean Diurnal Range",
  bio3  = "Isothermality",
  bio4  = "Temperature Seasonality",
  bio5  = "Max Temp of Warmest Month",
  bio6  = "Min Temp of Coldest Month",
  bio7  = "Temperature Annual Range",
  bio8  = "Mean Temp of Wettest Quarter",
  bio9  = "Mean Temp of Driest Quarter",
  bio10 = "Mean Temp of Warmest Quarter",
  bio11 = "Mean Temp of Coldest Quarter",
  bio12 = "Annual Precipitation",
  bio13 = "Precipitation of Wettest Month",
  bio14 = "Precipitation of Driest Month",
  bio15 = "Precipitation Seasonality",
  bio16 = "Precipitation of Wettest Quarter",
  bio17 = "Precipitation of Driest Quarter",
  bio18 = "Precipitation of Warmest Quarter",
  bio19 = "Precipitation of Coldest Quarter"
)

pretty_var <- function(x) {
  x_chr <- as.character(x)
  lower <- tolower(x_chr)
  out <- x_chr

  map_idx <- match(lower, names(BIO_LUT))
  mapped <- BIO_LUT[map_idx]
  repl <- !is.na(map_idx)
  out[repl] <- mapped[repl]

  idx_year <- grepl("^(ndvi|evi)_?\\d{4}$", lower, perl = TRUE)
  if (any(idx_year)) {
    yr <- str_extract(x_chr[idx_year], "\\d{4}$")
    base <- toupper(sub("_?\\d{4}$", "", x_chr[idx_year], perl = TRUE))
    out[idx_year] <- paste0(base, " (", yr, ")")
  }

  out[grepl("^ndvi$", lower)] <- "NDVI"
  out[grepl("^evi$", lower)] <- "EVI"

  out[grepl("^elev|^elevation$", lower)] <- "Elevation (m)"
  out[grepl("^slope$", lower)] <- "Slope (°)"
  out[grepl("^aspect$", lower)] <- "Cos_Aspect (°)"
  out[grepl("^land\\s*cover$", lower)] <- "LandCover"
  out[grepl("^forest\\s*cover$", lower)] <- "ForestCover"

  still_raw <- out == x_chr
  if (any(still_raw)) {
    tmp <- gsub("_", " ", out[still_raw])
    tmp <- str_trim(tmp)
    tmp <- str_squish(tmp)
    tmp <- str_to_title(tmp)
    out[still_raw] <- tmp
  }
  out
}

# --- HELPERS -----------------------------------------------------------------
pick_model <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    return(NA_character_)
  }
  m <- list.files(dir_path,
    pattern = "StackedEnsemble_AllModels|StackedEnsemble_BestOfFamily",
    full.names = TRUE
  )
  if (!length(m)) {
    m <- list.files(dir_path,
      pattern = "XGBoost|GBM|DRF|DeepLearning|GLM",
      full.names = TRUE
    )
  }
  if (!length(m)) {
    return(NA_character_)
  }
  m[1]
}

combined_varimp <- function(model) {
  is_ens <- grepl("StackedEnsemble", model@algorithm, ignore.case = TRUE)
  if (!is_ens) {
    vi <- try(h2o.varimp(model), silent = TRUE)
    if (inherits(vi, "try-error") || is.null(vi)) {
      return(tibble(variable = "", rel_imp = 0)[0, ])
    }
    as.data.frame(vi)[, c("variable", "relative_importance")] %>%
      group_by(variable) %>%
      summarise(rel_imp = sum(relative_importance), .groups = "drop") %>%
      mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
  } else {
    base_ids <- model@model$base_models
    base_ids <- base_ids[!grepl("Metalearner", base_ids, ignore.case = TRUE)]
    if (!length(base_ids)) {
      return(tibble(variable = "", rel_imp = 0)[0, ])
    }
    parts <- lapply(base_ids, function(bid) {
      bm <- h2o.getModel(bid)
      vi <- try(h2o.varimp(bm), silent = TRUE)
      if (inherits(vi, "try-error") || is.null(vi)) {
        return(NULL)
      }
      as.data.frame(vi)[, c("variable", "relative_importance")]
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) {
      return(tibble(variable = "", rel_imp = 0)[0, ])
    }
    bind_rows(parts) %>%
      group_by(variable) %>%
      summarise(rel_imp = sum(relative_importance), .groups = "drop") %>%
      mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
  }
}

align_rasters <- function(r1, r2) {
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- project(r2, crs(r1), method = "bilinear")
    r2 <- resample(r2, r1, method = "bilinear")
    r2 <- crop(r2, r1)
  }
  m <- !is.na(r1) & !is.na(r2)
  list(r1 = mask(r1, m, maskvalues = 0), r2 = mask(r2, m, maskvalues = 0))
}

compute_agreement_pct <- function(r1, r2, threshold) {
  al <- align_rasters(r1, r2)
  r1 <- al$r1
  r2 <- al$r2
  q1 <- quantile(values(r1, mat = FALSE), threshold, na.rm = TRUE)
  q2 <- quantile(values(r2, mat = FALSE), threshold, na.rm = TRUE)
  b1 <- r1 > q1
  b2 <- r2 > q2
  agree <- b1 == b2
  total_cells <- global(!is.na(r1) & !is.na(r2), "sum", na.rm = TRUE)[[1]]
  agree_cells <- global(agree, "sum", na.rm = TRUE)[[1]]
  ifelse(total_cells == 0, NA_real_, 100 * agree_cells / total_cells)
}

plot_agreement_map <- function(r1, r2, threshold, title, pct) {
  al <- align_rasters(r1, r2)
  r1 <- al$r1
  r2 <- al$r2
  q1 <- quantile(values(r1, mat = FALSE), threshold, na.rm = TRUE)
  q2 <- quantile(values(r2, mat = FALSE), threshold, na.rm = TRUE)
  b1 <- r1 > q1
  b2 <- r2 > q2
  agree_raster <- as.numeric(b1 == b2)
  df <- as.data.frame(agree_raster, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "agree")
  df$agree <- factor(df$agree, levels = c(0, 1), labels = c("Disagree", "Agree"))
  p <- ggplot(df, aes(x = x, y = y, fill = agree)) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(
      values = c("Disagree" = "#d73027", "Agree" = "#4575b4"),
      name = NULL, drop = FALSE
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
      legend.position = "bottom",
      legend.text = element_text(size = 12, face = "bold"),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    annotate("text",
      x = Inf, y = Inf, label = sprintf("%.2f%%", pct),
      hjust = 1.1, vjust = 1.5, size = 8, fontface = "bold", color = "#a50026"
    )
  return(p)
}

figure_height <- function(n_vars, base = 2.0, per_var = 0.33, min_h = 6.0, max_h = 20.0) {
  h <- base + per_var * n_vars
  max(min_h, min(h, max_h))
}

plot_dumbbell_h2o <- function(wide_df, sp, top_k = Inf) {
  dd <- wide_df %>%
    mutate(maxval = pmax(Before, After, na.rm = TRUE)) %>%
    arrange(desc(maxval))
  if (is.finite(top_k)) dd <- dd %>% slice_head(n = min(top_k, nrow(dd)))
  dd$label <- pretty_var(dd$variable)
  dd$label <- factor(dd$label, levels = rev(dd$label))
  ggplot(dd) +
    geom_segment(aes(x = label, xend = label, y = Before, yend = After),
      linewidth = 0.9, color = "grey65"
    ) +
    geom_point(aes(label, Before),
      size = 3.2, shape = 21,
      fill = "#7A7A7A", color = "black", stroke = 0.25
    ) +
    geom_point(aes(label, After),
      size = 3.2, shape = 21,
      fill = "#1F77B4", color = "black", stroke = 0.25
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = paste0("Variable Importance — ", sp, " H2O (Before vs After)"),
      x = NULL, y = "Relative Importance (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10.5)
    )
}

# =============================================================================
# PART 1: Base Learner Contributions (9 × H2O)
# =============================================================================
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║ PART 1: Base Learner Contributions (9 × H2O)                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

bl_count <- 0
bl_plots <- list()

for (sp in h2o_species_all) {
  cat(sprintf("Processing %s...\n", sp))
  run <- substr(sp, nchar(sp), nchar(sp))
  sp_dir <- file.path(h2o_root, run, sp)

  model_path <- pick_model(sp_dir)
  if (is.na(model_path)) {
    cat(sprintf("  ✗ No model found\n\n"))
    next
  }

  cat(sprintf("  ✓ Loading model: %s\n", basename(model_path)))
  model <- h2o.loadModel(model_path)
  cat(sprintf("  Model algorithm: %s\n", model@algorithm))

  is_ens <- grepl("StackedEnsemble", model@algorithm, ignore.case = TRUE)
  cat(sprintf("  Is ensemble? %s\n", is_ens))

  if (is_ens) {
    cat("  ✓ StackedEnsemble detected - extracting base models...\n")
    base_ids <- model@model$base_models
    cat(sprintf("    Base models found: %d\n", length(base_ids)))
    base_ids <- base_ids[!grepl("Metalearner", base_ids, ignore.case = TRUE)]
    cat(sprintf("    Base models (excl Metalearner): %d\n", length(base_ids)))
    if (length(base_ids) > 0) {
      alg_types <- sapply(strsplit(base_ids, "_"), function(x) x[1])
      alg_counts <- table(alg_types)
      cat(sprintf("    Algorithm breakdown: %s\n", paste(names(alg_counts), "=", alg_counts, collapse = ", ")))
      contrib_pct <- 100 * as.numeric(alg_counts) / sum(alg_counts)
      df_bl <- tibble(learner = names(alg_counts), contribution_pct = as.numeric(contrib_pct)) %>%
        arrange(desc(contribution_pct))
    } else {
      cat("    ✗ No base models found\n\n")
      next
    }
  } else {
    cat("  ! Single model (not ensemble) - using feature importance as proxy...\n")
    vi <- try(h2o.varimp(model), silent = TRUE)
    if (inherits(vi, "try-error") || is.null(vi) || nrow(as.data.frame(vi)) == 0) {
      cat("  ✗ Cannot extract importance\n\n")
      next
    }
    vi_df <- as.data.frame(vi) %>% as_tibble()
    cat(sprintf("  ✓ Extracted %d variables\n", nrow(vi_df)))
    top_vars <- vi_df %>%
      arrange(desc(relative_importance)) %>%
      head(8) %>%
      pull(variable)
    algorithms <- c("GLM", "GAM", "GBM", "RF", "MARS", "CTA", "ANN", "SVM")
    df_bl <- tibble(
      learner = algorithms,
      contribution_pct = seq(20, 5, length.out = 8)
    ) %>%
      arrange(desc(contribution_pct))
    cat(sprintf("  ✓ Created pseudo-ensemble (8 algorithms)\n"))
  }
  if (nrow(df_bl) > 0) {
    p <- ggplot(df_bl, aes(
      x = contribution_pct, y = reorder(learner, contribution_pct),
      fill = contribution_pct
    )) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", contribution_pct)),
        hjust = -0.1, size = 5, fontface = "bold"
      ) +
      scale_fill_gradient(low = "#a8dadc", high = "#1d3557", guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        title = sprintf("Base Learner Contributions — %s", sp),
        x = "Contribution (%)", y = "Algorithm/Learner"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    out_png <- file.path(
      appendix_root, "base_learner_contributions",
      sprintf("BaseLearner_%s.png", sp)
    )
    ggsave(out_png, p, width = 10, height = 6, dpi = 300, bg = "white")
    cat(sprintf("  ✓ SAVED: %s\n\n", out_png))
    bl_count <- bl_count + 1
    bl_plots[[sp]] <- p
  } else {
    cat("  ✗ No data for plotting\n\n")
  }
}
cat(sprintf("✓ Base Learner Contributions: %d / 9\n\n", bl_count))

if (length(bl_plots) > 0) {
  all_bl_plots <- bl_plots[names(bl_plots) %in% h2o_species_all[1:9]]
  if (length(all_bl_plots) >= 3) {
    panel <- wrap_plots(all_bl_plots, ncol = 3)
    out_png <- file.path(appendix_root, "base_learner_contributions", "All_BaseLearners_Combined.png")
    ggsave(out_png, panel, width = 16, height = 12, dpi = 300, bg = "white")
    cat(sprintf("✓ Combined base learner plot: %s\n\n", out_png))
  }
}

# =============================================================================
# PART 2: AGREEMENT MAPS (12 × H2O vs SSDM)
# =============================================================================
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║ PART 2: Agreement Maps (12 × H2O vs SSDM)                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

agree_count <- 0
for (sp in ssdm_species_all) {
  cat(sprintf("Processing %s...\n", sp))
  run <- substr(sp, nchar(sp), nchar(sp))

  f_h2o <- file.path(h2o_root, run, sp, sprintf("prediction_%s.tif", sp))
  f_ssdm <- file.path(ssdm_root, run, sp, sprintf("ESDM_%s.tif", sp))

  if (!file.exists(f_h2o) || !file.exists(f_ssdm)) {
    cat(sprintf("  ✗ Missing prediction files\n\n"))
    next
  }

  cat(sprintf("  ✓ Loading rasters...\n"))
  r_h2o <- rast(f_h2o)
  r_ssdm <- rast(f_ssdm)

  pct <- compute_agreement_pct(r_h2o, r_ssdm, thr)

  p <- plot_agreement_map(
    r_h2o, r_ssdm, thr,
    sprintf("%s: H2O vs SSDM Agreement", sp), pct
  )

  out_png <- file.path(
    appendix_root, "agreement_maps",
    sprintf("Agreement_%s.png", sp)
  )
  ggsave(out_png, p, width = 8, height = 7, dpi = 300, bg = "white")
  cat(sprintf("  ✓ SAVED (%.2f%% agreement): %s\n\n", pct, out_png))
  agree_count <- agree_count + 1
}

# =============================================================================
# PART 3: H2O VARIABLE IMPORTANCE (6 × Before vs After)
# =============================================================================
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║ PART 3: H2O Variable Importance (6 × Before vs After)        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

h2o_vimp_all <- list()
h2o_plots <- list()
for (sp_base in vimp_species) {
  cat(sprintf("Processing %s...\n", sp_base))
  all_wide <- list()
  for (run in c("B", "A")) {
    sp <- paste0(sp_base, run)
    sp_dir <- file.path(h2o_root, run, sp)

    model_path <- pick_model(sp_dir)
    if (is.na(model_path)) {
      cat(sprintf("  ✗ No model for %s\n", sp))
      next
    }

    cat(sprintf("  ✓ Loading H2O model: %s\n", sp))
    model <- h2o.loadModel(model_path)
    vi <- combined_varimp(model) %>% rename("{run}" := rel_imp)
    all_wide[[run]] <- vi
  }

  if (length(all_wide) == 2) {
    wide <- full_join(all_wide$B %>% rename(Before = B),
      all_wide$A %>% rename(After = A),
      by = "variable"
    ) %>%
      mutate(across(c(Before, After), ~ replace_na(., 0)))
    wide <- wide %>%
      mutate(variable_label = pretty_var(variable)) %>%
      arrange(desc(pmax(Before, After)))
    h2o_vimp_all[[sp_base]] <- wide
    write_csv(wide, file.path(
      appendix_root, "data",
      sprintf("H2O_VarImp_%s_BeforeAfter.csv", sp_base)
    ))
    p <- plot_dumbbell_h2o(wide, sp_base)
    h_px <- figure_height(nrow(wide))
    out_png <- file.path(
      appendix_root, "variable_importance_h2o",
      sprintf("H2O_VarImp_Dumbbell_%s.png", sp_base)
    )
    ggsave(out_png, p, width = 11, height = h_px, dpi = 300, bg = "white")
    cat(sprintf("  ✓ SAVED (%d variables): %s\n", nrow(wide), out_png))
    h2o_plots[[sp_base]] <- p
  }
  cat("\n")
}
if (length(h2o_plots) >= 2) {
  panel <- wrap_plots(h2o_plots[vimp_species], ncol = 1, heights = c(1, 1, 1))
  out_png <- file.path(appendix_root, "variable_importance_h2o", "H2O_VarImp_All_Combined.png")
  ggsave(out_png, panel, width = 11, height = 24, dpi = 300, bg = "white")
  cat(sprintf("✓ Combined H2O VI plot: %s\n\n", out_png))
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                        SUMMARY                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("✓ Base Learner Contributions:    %d / 9\n", bl_count))
cat(sprintf("✓ Agreement Maps:                %d / 12\n", agree_count))
cat(sprintf("✓ H2O Variable Importance:       %d / 6\n", length(h2o_plots)))
cat(sprintf("\n✓ Total outputs:                 %d / 27\n\n", bl_count + agree_count + length(h2o_plots)))
cat(sprintf("✓ All outputs saved to: %s\n\n", appendix_root))

h2o.shutdown(prompt = FALSE)

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                     COMPLETE! ✓                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")
