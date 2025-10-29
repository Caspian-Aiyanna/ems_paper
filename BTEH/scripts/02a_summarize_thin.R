#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fs)
  library(stringr)
  library(tidyr)
  library(grid)
  library(gridExtra)
  library(yaml)
  library(tibble)
  library(scales)
})

# --- paths (relative to project root) ---
dir_raw   <- "data/clean"
dir_main  <- "data/occ/thinned_DBSCAN"
dir_unc   <- "replicates_DBSCAN"              # <-- matches 02_dbscan_thin.R output
plans_yml <- "plans/dbscan_params.yml"

out_dir   <- "results/compare"
out_png   <- fs::path(out_dir, "tables.png")
out_pdf   <- fs::path(out_dir, "tables.pdf")
out_csv   <- fs::path(out_dir, "thin_summary.csv")

if (!dir_exists(out_dir)) dir_create(out_dir, recurse = TRUE)

# ---- Load run metadata (optional) ----
run_meta <- list(run = NA_character_, eps_deg = NA_real_, minPts = NA_integer_,
                 fraction = NA_real_, min_samples = NA_integer_, timestamp = NA_character_)
if (file_exists(plans_yml)) {
  prm <- tryCatch(yaml::read_yaml(plans_yml), error = function(e) NULL)
  if (!is.null(prm) && length(prm)) {
    pick <- tryCatch({
      runs <- names(prm)
      if (length(runs)) {
        ts <- sapply(runs, function(r) prm[[r]]$dbscan$timestamp %||% "")
        runs[order(ts, decreasing = TRUE)][1]
      } else NA_character_
    }, error = function(e) NA_character_)
    if (!is.na(pick)) {
      b <- prm[[pick]]$dbscan
      run_meta <- list(
        run        = pick,
        eps_deg    = b$eps_deg,
        minPts     = b$minPts,
        fraction   = b$fraction,
        min_samples= b$min_samples,
        timestamp  = b$timestamp
      )
    }
  }
}

# ---- helpers ----
species_from_path <- function(p) fs::path_ext_remove(fs::path_file(p))
safe_nrow <- function(p) {
  if (!file_exists(p)) return(NA_integer_)
  suppressMessages(suppressWarnings(nrow(readr::read_csv(p, show_col_types = FALSE))))
}
fmt_int  <- function(x) ifelse(is.na(x), "—", scales::comma(x, accuracy = 1))
fmt_pct1 <- function(x) ifelse(is.na(x), "—", sprintf("%.1f%%", x))
fmt_num1 <- function(x) ifelse(is.na(x), "—", sprintf("%.1f", x))

# ---- species list from MAIN outputs ----
main_files <- dir_ls(dir_main, glob = "*.csv", recurse = FALSE)
if (!length(main_files)) stop("No CSVs found in data/occ/thinned_DBSCAN — run 02_dbscan_thin.R first.")
species <- sort(species_from_path(main_files))

# ---- build main summary ----
main_df <- tibble(
  species = species,
  raw_n   = vapply(fs::path(dir_raw, paste0(species, ".csv")), safe_nrow, integer(1)),
  main_n  = vapply(fs::path(dir_main, paste0(species, ".csv")), safe_nrow, integer(1))
) |>
  mutate(keep_pct = if_else(is.finite(main_n / pmax(1, raw_n)),
                            100 * main_n / pmax(1, raw_n), NA_real_))

# ---- build uncertainty summary (rep1..rep3) ----
unc_df <- tibble(
  species = species,
  rep1_n  = vapply(fs::path(dir_unc, paste0(species, "_rep1.csv")), safe_nrow, integer(1)),
  rep2_n  = vapply(fs::path(dir_unc, paste0(species, "_rep2.csv")), safe_nrow, integer(1)),
  rep3_n  = vapply(fs::path(dir_unc, paste0(species, "_rep3.csv")), safe_nrow, integer(1))
) |>
  rowwise() |>
  mutate(
    rep_mean = mean(c_across(starts_with("rep")), na.rm = TRUE),
    rep_sd   = sd(  c_across(starts_with("rep")), na.rm = TRUE),
    rep_cv   = ifelse(is.finite(rep_mean) && rep_mean > 0, 100 * rep_sd / rep_mean, NA_real_)
  ) |>
  ungroup()

# ---- write CSV for downstream use ----
summary_csv <- main_df |>
  left_join(unc_df, by = "species") |>
  arrange(species)
readr::write_csv(summary_csv, out_csv)

# ---- pretty tables ----
main_print <- main_df |>
  transmute(
    Species   = species,
    `Raw (n)` = fmt_int(raw_n),
    `Main (n)`= fmt_int(main_n),
    `Keep %`  = fmt_pct1(keep_pct)
  )
tot_row1 <- tibble(
  Species   = "Total",
  `Raw (n)` = fmt_int(sum(main_df$raw_n, na.rm = TRUE)),
  `Main (n)`= fmt_int(sum(main_df$main_n, na.rm = TRUE)),
  `Keep %`  = fmt_pct1(100 * sum(main_df$main_n, na.rm = TRUE) / pmax(1, sum(main_df$raw_n, na.rm = TRUE)))
)
main_print <- bind_rows(main_print, tot_row1)

unc_print <- unc_df |>
  transmute(
    Species    = species,
    `Rep1 (n)` = fmt_int(rep1_n),
    `Rep2 (n)` = fmt_int(rep2_n),
    `Rep3 (n)` = fmt_int(rep3_n),
    `Mean (n)` = fmt_int(rep_mean),
    `SD`       = fmt_num1(rep_sd),
    `CV %`     = fmt_num1(rep_cv)
  )
tot_row2 <- tibble(
  Species    = "Total",
  `Rep1 (n)` = fmt_int(sum(unc_df$rep1_n, na.rm = TRUE)),
  `Rep2 (n)` = fmt_int(sum(unc_df$rep2_n, na.rm = TRUE)),
  `Rep3 (n)` = fmt_int(sum(unc_df$rep3_n, na.rm = TRUE)),
  `Mean (n)` = fmt_int(sum(unc_df$rep_mean, na.rm = TRUE)),
  `SD`       = fmt_num1(sqrt(sum(unc_df$rep_sd^2, na.rm = TRUE))),
  `CV %`     = "—"
)
unc_print <- bind_rows(unc_print, tot_row2)

tt <- ttheme_minimal(
  core = list(
    bg_params = list(fill = rep(c("#FFFFFF", "#F6F6F6"), length.out = nrow(main_print))),
    fg_params = list(cex = 0.9)
  ),
  colhead = list(fg_params = list(fontface = "bold", cex = 0.95))
)
tt2 <- ttheme_minimal(
  core = list(
    bg_params = list(fill = rep(c("#FFFFFF", "#F6F6F6"), length.out = nrow(unc_print))),
    fg_params = list(cex = 0.9)
  ),
  colhead = list(fg_params = list(fontface = "bold", cex = 0.95))
)

tab1 <- tableGrob(main_print, rows = NULL, theme = tt)
tab2 <- tableGrob(unc_print,  rows = NULL, theme = tt2)

title_text <- "DBSCAN Thinning Summary"
sub_text <- sprintf(
  "Run: %s | eps: %s deg | minPts: %s | fraction: %s | min_samples: %s | %s",
  run_meta$run %||% "—",
  ifelse(is.na(run_meta$eps_deg), "—", run_meta$eps_deg),
  ifelse(is.na(run_meta$minPts), "—", run_meta$minPts),
  ifelse(is.na(run_meta$fraction), "—", run_meta$fraction),
  ifelse(is.na(run_meta$min_samples), "—", run_meta$min_samples),
  run_meta$timestamp %||% format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
title_grob <- textGrob(title_text, gp = gpar(fontsize = 16, fontface = "bold"))
sub_grob   <- textGrob(sub_text, gp = gpar(fontsize = 9, col = "gray30"))
cap_grob   <- textGrob(
  "Notes: ‘Keep %’ = Main / Raw. CV %% is across replicate counts (higher = less stable). ‘—’ denotes not available.",
  gp = gpar(fontsize = 8, col = "gray35")
)
h1 <- textGrob("Main replicate", gp = gpar(fontface = "bold"))
h2 <- textGrob("Uncertainty replicates", gp = gpar(fontface = "bold"))

g <- arrangeGrob(
  title_grob, sub_grob,
  arrangeGrob(h1, tab1, ncol = 1, heights = c(unit(1, "lines"), unit(1, "null"))),
  arrangeGrob(h2, tab2, ncol = 1, heights = c(unit(1, "lines"), unit(1, "null"))),
  cap_grob,
  ncol = 1,
  heights = unit.c(unit(1.2, "lines"), unit(0.9, "lines"), unit(1, "null"), unit(1, "null"), unit(0.8, "lines"))
)

# --- Save PNG (300 dpi) ---
png(out_png, width = 2000, height = 2600, res = 300)
grid.draw(g)
dev.off()

# --- Save PDF (vector) ---
pdf(out_pdf, width = 8.5, height = 11)
grid.draw(g)
dev.off()

cat("Wrote:\n -", out_png, "\n -", out_pdf, "\n -", out_csv, "\n")
