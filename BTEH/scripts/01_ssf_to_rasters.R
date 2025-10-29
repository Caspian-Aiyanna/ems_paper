suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(amt)
  library(sf)
  library(terra)
  library(survival)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

# ============================ USER CONFIG =====================================
project_root <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/BTEH/BTEH"

dirs <- list(
  A = list(env_dir = "data/envi/A/stack", sp_dir = "data/clean/A", res_dir = "results/SSF/A"),
  B = list(env_dir = "data/envi/B/stack", sp_dir = "data/clean/B", res_dir = "results/SSF/B")
)

tz_data     <- "UTC"
seed_val    <- 20161113
n_controls  <- 10             # random steps per used step
min_burst_n <- 4              # ≥(min_burst_n) fixes per burst
min_steps_for_fit <- 20       # minimum used steps to estimate Gamma
drop_na_prop <- 0.3           # drop covariates with >30% NA at endpoints
near_zero_sd <- 1e-8          # drop near-constant covariates

# Plots
max_vars_partial <- 6         # number of top-|beta| variables for partial response plots
fig_width <- 9; fig_height <- 6; dpi <- 300

set.seed(seed_val)

# ============================ HELPERS =========================================
req_cols <- c("collar_id", "species", "timestamp", "lon", "lat")

parse_ts <- function(x, tz = "UTC") {
  ts <- suppressWarnings(lubridate::dmy_hm(x, tz = tz))
  if (all(is.na(ts))) ts <- suppressWarnings(lubridate::ymd_hms(x, tz = tz))
  if (all(is.na(ts))) ts <- suppressWarnings(lubridate::ymd_hm(x, tz = tz))
  if (all(is.na(ts))) stop("Timestamp parsing failed. Check 'timestamp' format.")
  ts
}

# Generate control endpoints (uniform angle; gamma step length)
make_controls <- function(steps_used, n_controls, gamma_shape, gamma_scale) {
  need <- c("x1_", "y1_", "step_id_")
  if (!all(need %in% names(steps_used))) {
    stop("Expected columns x1_, y1_, step_id_ in steps table. Have: ",
         paste(names(steps_used), collapse = ", "))
  }
  N  <- nrow(steps_used)
  if (N == 0L) return(dplyr::slice_head(steps_used, n = 0))

  base <- steps_used[rep(seq_len(N), each = n_controls), c("x1_", "y1_", "step_id_")]
  L    <- rgamma(n = N * n_controls, shape = gamma_shape, scale = gamma_scale)
  ang  <- runif(n = N * n_controls, min = 0, max = 2*pi)

  tibble(
    step_id_ = base$step_id_,
    x2_ = base$x1_ + L * cos(ang),
    y2_ = base$y1_ + L * sin(ang),
    case_ = 0L
  )
}

# Dynamic band names one..twentysix then band27+
num2word <- function(k){
  w <- c("one","two","three","four","five","six","seven","eight","nine","ten",
         "eleven","twelve","thirteen","fourteen","fifteen","sixteen","seventeen",
         "eighteen","nineteen","twenty","twentyone","twentytwo","twentythree",
         "twentyfour","twentyfive","twentysix")
  if (k <= length(w)) w[k] else paste0("band", k)
}

pretty_var <- function(v) gsub("_", " ", v)

theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# ============================ CORE PER-FILE PROCESS ===========================
process_one_file <- function(run_tag, stack_path, csv_file, out_base, i_seed = 0L) {
  tag <- tools::file_path_sans_ext(basename(csv_file))
  message("\n--- [", run_tag, "] Processing: ", tag, " ---")

  out_dir  <- file.path(out_base, tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "run_log.txt")

  log_msg <- function(...) {
    txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    message(txt); cat(txt, "\n", file = log_file, append = TRUE)
  }

  tryCatch({
    # ---------- Load & prepare stack ----------
    stopifnot(file.exists(stack_path))
    r_stack    <- terra::rast(stack_path)
    stack_crs  <- sf::st_crs(r_stack)

    # Save original names, then assign safe sequential names
    orig_names <- names(r_stack)
    n_bands <- nlyr(r_stack)
    new_names <- vapply(seq_len(n_bands), num2word, character(1))
    names(r_stack) <- make.names(new_names, unique = TRUE)
    safe_names <- names(r_stack)

    # Save name mapping (original -> safe) for traceability
    readr::write_csv(
      tibble(original = orig_names, safe = safe_names),
      file.path(out_dir, "stack_name_mapping.csv")
    )

    # ---------- Load data ----------
    stopifnot(file.exists(csv_file))
    dat <- readr::read_csv(csv_file, show_col_types = FALSE)

    miss <- setdiff(req_cols, names(dat))
    if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

    dat$timestamp <- parse_ts(dat$timestamp, tz = tz_data)
    if (any(is.na(dat$lon) | is.na(dat$lat))) stop("Missing lon/lat values.")
    if (any(abs(dat$lon) > 180 | abs(dat$lat) > 90)) stop("Invalid lon/lat values.")

    dat <- dat %>%
      arrange(collar_id, timestamp) %>%
      group_by(collar_id) %>%
      filter(!(duplicated(timestamp) & duplicated(lon) & duplicated(lat))) %>%
      ungroup()

    pts_wgs <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)
    pts_m   <- sf::st_transform(pts_wgs, crs = stack_crs)
    crd     <- sf::st_coordinates(pts_m)
    pts_m$X <- crd[,1]; pts_m$Y <- crd[,2]

    trk <- amt::make_track(pts_m, .x = X, .y = Y, .t = timestamp, id = collar_id)

    # ---------- resample 30±10 -> fallback 60±20 ----------
    trk30 <- trk %>%
      amt::track_resample(rate = lubridate::minutes(30), tolerance = lubridate::minutes(10)) %>%
      amt::filter_min_n_burst(min_n = min_burst_n)

    trk_use <- if (nrow(trk30) >= 2 * min_burst_n) {
      log_msg("Using cadence 30±10 minutes.")
      trk30
    } else {
      log_msg("Insufficient fixes for 30±10; trying 60±20.")
      trk60 <- trk %>%
        amt::track_resample(rate = lubridate::minutes(60), tolerance = lubridate::minutes(20)) %>%
        amt::filter_min_n_burst(min_n = min_burst_n)
      if (nrow(trk60) < 2 * min_burst_n) stop("Insufficient fixes for both cadences; skipping.")
      log_msg("Using cadence 60±20 minutes.")
      trk60
    }

    # ---------- used steps ----------
    steps_raw <- amt::steps_by_burst(trk_use)

    if (!"step_id_" %in% names(steps_raw)) {
      if ("burst_" %in% names(steps_raw)) {
        steps_raw <- steps_raw %>%
          group_by(burst_) %>%
          mutate(step_id_ = paste0(as.character(burst_), "_", row_number())) %>%
          ungroup()
      } else {
        steps_raw$step_id_ <- sprintf("s%07d", seq_len(nrow(steps_raw)))
      }
      log_msg("Added missing step_id_ column.")
    }

    eps <- 1e-5
    if ("sl_" %in% names(steps_raw)) {
      zero_n <- sum(is.finite(steps_raw$sl_) & steps_raw$sl_ == 0)
      if (zero_n > 0) {
        log_msg("Found ", zero_n, " zero-length steps -> replacing with epsilon.")
        steps_raw$sl_[is.finite(steps_raw$sl_) & steps_raw$sl_ == 0] <- eps
      }
    }

    needed_cols <- c("x1_", "y1_", "x2_", "y2_")
    miss_cols <- setdiff(needed_cols, names(steps_raw))
    if (length(miss_cols) > 0) stop("Expected columns missing in steps_by_burst output: ", paste(miss_cols, collapse = ", "))

    # ---------- Gamma for controls ----------
    sl_vec <- steps_raw$sl_[is.finite(steps_raw$sl_) & steps_raw$sl_ > 0]
    if (length(sl_vec) < min_steps_for_fit) stop("Too few valid step lengths for param estimation.")
    sl_mean <- mean(sl_vec); sl_sd <- stats::sd(sl_vec)
    if (!is.finite(sl_sd) || sl_sd == 0) {
      log_msg("Zero SD for step lengths; inflating SD slightly.")
      sl_sd <- max(eps, sl_mean * 0.05)
    }
    gamma_shape <- (sl_mean / sl_sd)^2
    gamma_scale <- (sl_sd^2) / sl_mean

    used <- steps_raw %>% select(step_id_, x2_, y2_) %>% mutate(case_ = 1L)

    set.seed(seed_val + i_seed)
    controls <- make_controls(
      steps_used   = steps_raw,
      n_controls   = n_controls,
      gamma_shape  = gamma_shape,
      gamma_scale  = gamma_scale
    )

    ssf_data <- bind_rows(used, controls) %>%
      filter(is.finite(x2_), is.finite(y2_))
    if (nrow(ssf_data) < 50) stop("Too few steps for modeling (n<50).")

    # ---------- Extract ALL covariates at endpoints ----------
    xy_mat <- as.matrix(ssf_data[, c("x2_", "y2_")])
    ext    <- terra::extract(r_stack, xy_mat, method = "bilinear")
    if ("ID" %in% names(ext)) ext <- ext[, -1, drop = FALSE]
    colnames(ext) <- names(r_stack)
    covar_names <- colnames(ext)
    ssf_data <- bind_cols(ssf_data, as_tibble(ext))

    # ---------- Filter covariates by NA and variance ----------
    na_rate <- sapply(ssf_data[, covar_names, drop = FALSE], function(z) mean(!is.finite(z)))
    keep1   <- names(na_rate)[na_rate <= drop_na_prop]
    sds     <- sapply(ssf_data[, keep1, drop = FALSE], stats::sd, na.rm = TRUE)
    keep2   <- names(sds)[is.finite(sds) & sds > near_zero_sd]
    cov_keep <- intersect(keep1, keep2)
    if (length(cov_keep) < 2) stop("Too few usable covariates after NA/variance filtering.")

    cc <- stats::complete.cases(ssf_data[, cov_keep, drop = FALSE])
    ssf_data <- ssf_data[cc, , drop = FALSE]
    if (nrow(ssf_data) < 50) stop("Too few steps after covariate filtering.")

    # ---------- Z-score ----------
    stats_df <- tibble(
      variable = cov_keep,
      mean = sapply(ssf_data[, cov_keep, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
      sd   = sapply(ssf_data[, cov_keep, drop = FALSE], function(z) stats::sd(z, na.rm = TRUE))
    )
    stats_df$sd[!is.finite(stats_df$sd) | stats_df$sd == 0] <- 1

    for (v in stats_df$variable) {
      zname <- paste0(v, "_z")
      mu <- stats_df$mean[stats_df$variable == v]
      sdv <- stats_df$sd[stats_df$variable == v]
      ssf_data[[zname]] <- (ssf_data[[v]] - mu) / sdv
    }

    # ---------- Model ----------
    z_vars <- paste0(cov_keep, "_z")
    base_fml <- stats::reformulate(termlabels = z_vars, response = "case_")
    fml <- update(base_fml, . ~ . + strata(step_id_))
    model_clogit <- survival::clogit(fml, data = ssf_data)

    # ====== Hardened coefficient/CI assembly (name-aligned; no out-of-bounds) ======
    sm <- summary(model_clogit)
    betas <- coef(model_clogit)
    nm <- names(betas)

    # Coef table as data.frame, row-aligned to nm
    coef_mat <- if (is.null(dim(sm$coefficients))) {
      cm <- as.data.frame(t(sm$coefficients))
      rownames(cm) <- nm
      cm
    } else {
      as.data.frame(sm$coefficients)
    }
    # Ensure rows are in the same order as nm
    if (!is.null(rownames(coef_mat))) {
      coef_mat <- coef_mat[nm, , drop = FALSE]
    } else {
      rownames(coef_mat) <- nm
    }

    # Confidence intervals aligned to nm
    ci_raw <- suppressWarnings(confint(model_clogit))
    if (is.null(dim(ci_raw))) {
      ci_mat <- matrix(ci_raw, nrow = 1, dimnames = list(nm, c("2.5 %","97.5 %")))
    } else {
      ci_mat <- as.matrix(ci_raw)
      if (!is.null(rownames(ci_mat))) {
        ci_mat <- ci_mat[nm, , drop = FALSE]
      } else {
        rownames(ci_mat) <- nm
      }
    }

    se   <- if ("SE(coef)" %in% colnames(coef_mat)) coef_mat[ , "SE(coef)"] else rep(NA_real_, length(nm))
    zval <- if ("z"        %in% colnames(coef_mat)) coef_mat[ , "z"]        else rep(NA_real_, length(nm))
    pval <- if ("Pr(>|z|)" %in% colnames(coef_mat)) coef_mat[ , "Pr(>|z|)"] else rep(NA_real_, length(nm))

    beta_tbl <- tibble(
      variable = sub("_z$", "", nm),
      beta     = as.numeric(betas),
      se       = as.numeric(se),
      z        = as.numeric(zval),
      p        = as.numeric(pval),
      ci_lo    = as.numeric(ci_mat[,1]),
      ci_hi    = as.numeric(ci_mat[,2])
    ) %>%
      left_join(stats_df, by = "variable") %>%
      select(variable, beta, se, z, p, ci_lo, ci_hi, mean, sd) %>%
      arrange(desc(abs(beta)))

    readr::write_csv(beta_tbl, file.path(out_dir, "model_coefficients_and_scaling.csv"))
    capture.output(sm, file = file.path(out_dir, "clogit_summary.txt"))

    # ---------- RSF projection (guard missing bands) ----------
    avail <- intersect(beta_tbl$variable, names(r_stack))
    if (length(avail) == 0) stop("No model covariates found in raster stack for projection.")
    beta_tbl <- beta_tbl %>% filter(variable %in% avail)
    sel <- r_stack[[avail]]

    mu_vec <- beta_tbl$mean
    sd_vec <- beta_tbl$sd; sd_vec[!is.finite(sd_vec) | sd_vec == 0] <- 1
    b_vec  <- beta_tbl$beta

    rsf <- terra::app(sel, fun = function(v) {
      z <- (v - mu_vec) / sd_vec
      exp(sum(z * b_vec))
    })

    g <- terra::global(rsf, fun = c("min","max"), na.rm = TRUE)
    rsf01 <- if (isTRUE(all.equal(g[1,1], g[1,2]))) rsf * 0 else (rsf - g[1,1]) / (g[1,2] - g[1,1])

    out_tif   <- file.path(out_dir, paste0(tag, "_SSF_rsf.tif"))
    out_tif01 <- file.path(out_dir, paste0(tag, "_SSF_rsf_0to1.tif"))
    terra::writeRaster(rsf,   out_tif,   overwrite = TRUE, wopt = list(datatype="FLT4S", gdal="COMPRESS=LZW"))
    terra::writeRaster(rsf01, out_tif01, overwrite = TRUE, wopt = list(datatype="FLT4S", gdal="COMPRESS=LZW"))

    # ===================== DIAGNOSTIC PLOTS (journal-ready) ===================
    if (nrow(beta_tbl) > 0) {
      coef_df <- beta_tbl %>%
        mutate(variable_label = pretty_var(variable)) %>%
        arrange(beta) %>%
        mutate(variable_label = factor(variable_label, levels = variable_label))

      p_coef <- ggplot(coef_df, aes(x = beta, y = variable_label)) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0) +
        geom_point(size = 2) +
        labs(title = paste0("SSF coefficients (", run_tag, " — ", tag, ")"),
             x = "Coefficient (log-relative selection)", y = NULL) +
        theme_pub()
      ggsave(file.path(out_dir, "coef_forest.png"), p_coef, width = fig_width, height = fig_height, dpi = dpi)

      p_vi <- ggplot(coef_df, aes(x = reorder(variable_label, abs(beta)), y = abs(beta))) +
        geom_col() + coord_flip() +
        labs(title = paste0("Variable importance (|β|) — ", run_tag, " — ", tag),
             x = NULL, y = "|β|") +
        theme_pub()
      ggsave(file.path(out_dir, "var_importance.png"), p_vi, width = fig_width, height = fig_height, dpi = dpi)

      top_vars <- head(coef_df$variable, min(10, nrow(coef_df)))
      dens_df <- ssf_data %>%
        select(case_, all_of(paste0(top_vars, "_z"))) %>%
        mutate(case_ = factor(case_, levels = c(0,1), labels = c("available","used"))) %>%
        pivot_longer(cols = -case_, names_to = "variable", values_to = "value") %>%
        mutate(variable = sub("_z$", "", variable),
               variable = factor(pretty_var(variable), levels = pretty_var(top_vars)))

      p_dens <- ggplot(dens_df, aes(x = value, fill = case_)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~ variable, scales = "free", ncol = 3) +
        labs(title = paste0("Endpoint covariate densities (z) — ", run_tag, " — ", tag),
             x = "z-score", y = "Density", fill = NULL) +
        theme_pub()
      ggsave(file.path(out_dir, "covariate_densities.png"), p_dens, width = fig_width, height = fig_height, dpi = dpi)

      pr_grid <- seq(-3, 3, length.out = 100)
      top_pr_vars <- head(coef_df$variable, min(max_vars_partial, nrow(coef_df)))
      pr_list <- lapply(top_pr_vars, function(v) {
        lp <- pr_grid * coef_df$beta[coef_df$variable == v]
        tibble(variable = pretty_var(v), z = pr_grid, RSF_rel = exp(lp))
      })
      pr_df <- bind_rows(pr_list)

      p_pr <- ggplot(pr_df, aes(x = z, y = RSF_rel)) +
        geom_line() +
        facet_wrap(~ variable, scales = "free_y", ncol = 3) +
        labs(title = paste0("Partial responses (others at mean) — ", run_tag, " — ", tag),
             x = "z-score of focal variable", y = "RSF (relative)") +
        theme_pub()
      ggsave(file.path(out_dir, "partial_responses.png"), p_pr, width = fig_width, height = fig_height, dpi = dpi)

      # Win rate
      Xz <- as.matrix(ssf_data[, paste0(coef_df$variable, "_z"), drop = FALSE])
      lp <- drop(Xz %*% matrix(coef_df$beta, ncol = 1))
      pred_rsf <- exp(lp)

      wr_df <- tibble(step_id_ = ssf_data$step_id_, case_ = ssf_data$case_, rsf = pred_rsf) %>%
        group_by(step_id_) %>%
        summarise(
          used_rsf = rsf[case_ == 1][1],
          max_rsf = max(rsf),
          win = as.integer(used_rsf >= max_rsf),
          .groups = "drop"
        )
      win_rate <- mean(wr_df$win, na.rm = TRUE)
      readr::write_csv(wr_df, file.path(out_dir, "winrate_per_step.csv"))

      p_wr <- ggplot(wr_df, aes(x = win)) +
        geom_bar() +
        scale_x_continuous(breaks = c(0,1), labels = c("lost","won")) +
        labs(title = paste0("Within-stratum win rate = ", percent(win_rate), " — ", run_tag, " — ", tag),
             x = NULL, y = "Count of strata") +
        theme_pub()
      ggsave(file.path(out_dir, "winrate_bar.png"), p_wr, width = 6, height = 4, dpi = dpi)
    }

    # ---------- Meta ----------
    meta <- tibble(
      run = run_tag,
      tag, input_csv = csv_file, n_points = nrow(dat),
      n_steps_used = sum(ssf_data$case_ == 1L),
      n_controls = n_controls,
      cadence = if (exists("trk30") && identical(trk_use, trk30)) "30±10" else "60±20",
      timestamp_tz = tz_data, seed = seed_val + i_seed,
      raster_stack = stack_path, output_rsf = out_tif, output_rsf01 = out_tif01,
      n_covariates = nrow(beta_tbl),
      covariates = paste(beta_tbl$variable, collapse = ", "),
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    if (exists("wr_df")) meta$win_rate <- mean(wr_df$win, na.rm = TRUE)

    readr::write_csv(meta, file.path(out_dir, "run_meta.csv"))
    log_msg("SUCCESS → ", out_tif, " (+ 0–1, plots, meta)")

  }, error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n", file = log_file, append = TRUE)
    message("ERROR in [", run_tag, "] ", tag, ": ", conditionMessage(e))
  })
}

# ============================ MAIN LOOP (A & B) ===============================
message("Project root: ", project_root)
for (k in names(dirs)) {
  test_stack <- file.path(project_root, dirs[[k]]$env_dir, "stack.tif")
  message("[", k, "] Stack candidate: ", test_stack, "  exists=", file.exists(test_stack))
}

for (run_tag in names(dirs)) {
  env_dir <- file.path(project_root, dirs[[run_tag]]$env_dir)
  sp_dir  <- file.path(project_root, dirs[[run_tag]]$sp_dir)
  out_dir_base <- file.path(project_root, dirs[[run_tag]]$res_dir)
  dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

  stack_path <- file.path(env_dir, "stack.tif")
  if (!file.exists(stack_path)) {
    message("[", run_tag, "] Missing stack: ", stack_path, " — skipping run.")
    next
  }

  csv_files <- list.files(sp_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    message("[", run_tag, "] No CSV files in: ", sp_dir, " — skipping run.")
    next
  }

  for (i in seq_along(csv_files)) {
    process_one_file(run_tag, stack_path, csv_files[[i]], out_dir_base, i_seed = i)
  }
}

message("\nAll done. Outputs under:\n  A → ", file.path(project_root, dirs$A$res_dir),
        "\n  B → ", file.path(project_root, dirs$B$res_dir))
