#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fs)
  library(stringr)
  library(tidyr)
  library(grid)
  library(gridExtra)
})

# --- paths (relative to project root)
dir_raw   <- "data/clean"
dir_main  <- "data/occ/thinned_DBSCAN"
dir_unc   <- "data/occ/uncertainty"
out_dir   <- "results/compare"
out_png   <- fs::path(out_dir, "tables.png")
out_csv   <- fs::path(out_dir, "thin_summary.csv")

if (!dir_exists(out_dir)) dir_create(out_dir, recurse = TRUE)

# helpers
species_from_path <- function(p) fs::path_ext_remove(fs::path_file(p))
safe_nrow <- function(p) {
  if (!file_exists(p)) return(NA_integer_)
  suppressMessages(suppressWarnings(nrow(readr::read_csv(p, show_col_types = FALSE))))
}

# species list from main outputs
main_files <- dir_ls(dir_main, glob = "*.csv", recurse = FALSE)
if (!length(main_files)) {
  stop("No CSVs found in data/occ/thinned_DBSCAN — run 02_dbscan_thin.R first.")
}
species <- sort(species_from_path(main_files))

# build main summary
main_df <- tibble(
  species = species,
  raw_n   = vapply(fs::path(dir_raw, paste0(species, ".csv")), safe_nrow, integer(1)),
  main_n  = vapply(fs::path(dir_main, paste0(species, ".csv")), safe_nrow, integer(1))
) |>
  mutate(keep_pct = if_else(is.finite(main_n / pmax(1, raw_n)),
                            round(100 * main_n / pmax(1, raw_n), 1), NA_real_))

# build uncertainty summary (rep1..rep3)
unc_df <- tibble(
  species = species,
  rep1_n  = vapply(fs::path(dir_unc, paste0(species, "_rep1.csv")), safe_nrow, integer(1)),
  rep2_n  = vapply(fs::path(dir_unc, paste0(species, "_rep2.csv")), safe_nrow, integer(1)),
  rep3_n  = vapply(fs::path(dir_unc, paste0(species, "_rep3.csv")), safe_nrow, integer(1))
)

# write a CSV for downstream use
summary_csv <- main_df |>
  left_join(unc_df, by = "species") |>
  arrange(species)
readr::write_csv(summary_csv, out_csv)

# pretty labels
main_print <- main_df |>
  transmute(
    Species = species,
    `Raw (n)` = raw_n,
    `Main (n)` = main_n,
    `Keep %` = sprintf("%.1f%%", keep_pct)
  )

unc_print <- unc_df |>
  transmute(
    Species = species,
    `Rep1 (n)` = rep1_n,
    `Rep2 (n)` = rep2_n,
    `Rep3 (n)` = rep3_n
  )

# make grobs
title_grob <- textGrob("DBSCAN Thinning Summary", gp = gpar(fontsize = 16, fontface = "bold"))
sub_grob   <- textGrob("Main replicate & uncertainty replicates", gp = gpar(fontsize = 10, col = "gray30"))

tab1 <- tableGrob(main_print, rows = NULL, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.9)),
                                                                  colhead = list(fg_params = list(fontface = "bold"))))
tab2 <- tableGrob(unc_print,  rows = NULL, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.9)),
                                                                  colhead = list(fg_params = list(fontface = "bold"))))

# arrange and save as PNG
g <- arrangeGrob(
  title_grob,
  sub_grob,
  arrangeGrob(textGrob("Main replicate", gp = gpar(fontface = "bold")), tab1, ncol = 1, heights = c(1, 10)),
  arrangeGrob(textGrob("Uncertainty replicates", gp = gpar(fontface = "bold")), tab2, ncol = 1, heights = c(1, 10)),
  ncol = 1,
  heights = unit.c(unit(1, "lines"), unit(1, "lines"), unit(1, "null"), unit(1, "null"))
)

# choose a width that fits table nicely (auto height)
png(out_png, width = 1400, height = 900, res = 150)
grid.draw(g)
dev.off()

cat("Wrote:\n -", out_png, "\n -", out_csv, "\n")
