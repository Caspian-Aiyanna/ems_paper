#!/usr/bin/env Rscript
# =============================================================================
# 06a_uncertainty.R
# Uncertainty Decomposition (SSF–RSF vs H2O & SSDM), per elephant & run (A/B)
#
# Inputs (supported patterns; script tries them in order and uses the first hit)
#   H2O  → results/H2O/<A|B>/<E#A|E#B>/prediction_thinned_<E#A|E#B>.tif
#        → results/H2O/<A|B>/prediction_thinned_<E#A|E#B>.tif
#   SSDM → results/SSDM/<A|B>/<E#A|E#B>/ESDM_<E#A|E#B>.tif
#        → results/SSDM/<A|B>/ESDM_<E#A|E#B>.tif
#   SSF  → results/SSF/<A|B>/<E#A|E#B>/<E#A|E#B>_SSF_rsf_0to1.tif
#        → results/SSF/<A|B>/<E#A|E#B>_SSF_rsf_0to1.tif
#        → results/SSF/<E#A|E#B>_SSF_rsf_0to1.tif
#
# Outputs:
#   paper_results/uncertainty/<ELEPHANT>/{H2O,SSDM,SSF,combined,figures,tables}/...
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(viridisLite)
  library(rnaturalearth)
})

# ------------------------------- CLI -----------------------------------------
option_list <- list(
  make_option(c("-e","--elephant"), type="character", default=NULL,
              help="Elephants: E3 or E3,E4,E5 or ALL. (default set used if source() without args)"),
  make_option(c("-a","--tagA"), type="character", default="A", help="Before tag [default %default]"),
  make_option(c("-b","--tagB"), type="character", default="B", help="After tag [default %default]"),
  make_option(c("-r","--root"), type="character", default="results", help="Results root [default %default]"),
  make_option(c("-q","--qmask"), type="double", default=0.75, help="Quantile for hotspot masks [default %default]"),
  make_option(c("-o","--outdir"), type="character", default="paper_results/uncertainty", help="Output root [default %default]"),
  make_option(c("-s","--seed"), type="integer", default=20161113L, help="Seed [default %default]"),
  make_option(c("--skip_figs"), action="store_true", default=FALSE, help="Skip PNG rendering")
)
opt <- parse_args(OptionParser(option_list=option_list))
# ---- Coastline handling (non-interactive & optional) ----
# Default: try to draw coastlines if the data pkgs are present.
opt$coast <- TRUE

# Allow you to force-disable when sourcing by defining `opt_disable_coast <- TRUE` in your console.
if (exists("opt_disable_coast", inherits = TRUE) && isTRUE(opt_disable_coast)) {
  opt$coast <- FALSE
}

# If required pkgs are missing in the *active* renv library, turn coastlines off silently.
if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
    (!requireNamespace("rnaturalearthdata", quietly = TRUE) &&
     !requireNamespace("rnaturalearthhires", quietly = TRUE))) {
  opt$coast <- FALSE
}

# Honor absolute root when sourcing()
if (exists("opt_root", inherits = TRUE) && nzchar(opt_root)) {
  opt$root <- opt_root
  message("[ROOT] Using opt_root = ", normalizePath(opt$root, mustWork = FALSE))
}

# Friendly fallback when sourced() without args
if (is.null(opt$elephant)) {
  message("No --elephant provided; assuming interactive run via source().")
  opt$elephant <- "E3,E4,E5,E1,E2,E6"   # base IDs only
}

set.seed(opt$seed)
runA <- toupper(opt$tagA); runB <- toupper(opt$tagB)

# --------------------------- Compatibility builders ---------------------------
first_existing <- function(...) {
  cands <- c(...)
  hit <- which(file.exists(cands))[1]
  if (is.na(hit)) return(cands[1])  # fall back to first (useful for debugging prints)
  cands[hit]
}

# H2O: try with per-elephant subfolder first, then flat in A/B
p_h2o <- function(ele, run) {
  c1 <- file.path(opt$root, "H2O", run, sprintf("%s%s", ele, run),
                  sprintf("prediction_thinned_%s%s.tif", ele, run))
  c2 <- file.path(opt$root, "H2O", run,
                  sprintf("prediction_thinned_%s%s.tif", ele, run))
  first_existing(c1, c2)
}

# SSDM: with per-elephant subfolder first, then flat in A/B
p_ssdm <- function(ele, run) {
  c1 <- file.path(opt$root, "SSDM", run, sprintf("%s%s", ele, run),
                  sprintf("ESDM_%s%s.tif", ele, run))
  c2 <- file.path(opt$root, "SSDM", run,
                  sprintf("ESDM_%s%s.tif", ele, run))
  first_existing(c1, c2)
}

# SSF: support nested per-elephant folder, flat in A/B, or fully flat
p_ssf <- function(ele, run) {
  c1 <- file.path(opt$root, "SSF", run, sprintf("%s%s", ele, run),
                  sprintf("%s%s_SSF_rsf_0to1.tif", ele, run))
  c2 <- file.path(opt$root, "SSF", run,
                  sprintf("%s%s_SSF_rsf_0to1.tif", ele, run))
  c3 <- file.path(opt$root, "SSF",
                  sprintf("%s%s_SSF_rsf_0to1.tif", ele, run))
  first_existing(c1, c2, c3)
}

exists_all <- function(ele, run) {
  all(file.exists(p_h2o(ele, run)),
      file.exists(p_ssdm(ele, run)),
      file.exists(p_ssf(ele, run)))
}

# --------------------------- Output dir helper --------------------------------
mk_outdirs <- function(ele) {
  out_root <- file.path(opt$outdir, ele)
  dirs <- list(
    H2O=file.path(out_root,"H2O"), SSDM=file.path(out_root,"SSDM"),
    SSF=file.path(out_root,"SSF"), combined=file.path(out_root,"combined"),
    figures=file.path(out_root,"figures"), tables=file.path(out_root,"tables")
  )
  invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))
  dirs
}

# --------------------------- Coastline + plotting ----------------------------
.coast_cache <- new.env(parent = emptyenv())

get_coastline <- function(template_rast){
  # Skip completely if disabled or if packages/data not available.
  if (!isTRUE(opt$coast)) return(NULL)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(NULL)
  if (!requireNamespace("rnaturalearthdata", quietly = TRUE) &&
      !requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    return(NULL)
  }

  key <- paste0("crs:", as.character(terra::crs(template_rast)))
  if (!exists(key, envir = .coast_cache)) {
    coast_sf <- tryCatch(
      rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf"),
      error = function(e) NULL
    )
    if (is.null(coast_sf)) return(NULL)
    coast_v <- terra::vect(coast_sf)
    if (!is.na(terra::crs(template_rast))) {
      coast_v <- tryCatch(terra::project(coast_v, template_rast), error = function(e) coast_v)
    }
    assign(key, coast_v, envir = .coast_cache)
  }
  get(key, envir = .coast_cache)
}


save_raster_map <- function(r, title, fname,
                            palette=c("sequential","diverging","binary","tri"),
                            width_px=1600, height_px=1200, dpi=240){
  palette <- match.arg(palette)
  dir.create(dirname(fname), TRUE, TRUE)
  v <- values(r); v <- v[is.finite(v)]
  png(fname, width=width_px, height=height_px, res=dpi, bg="white")
  op <- par(no.readonly=TRUE); on.exit({par(op); dev.off()}, add=TRUE)
  par(mar=c(3.2,3.6,3.2,6))
  coast <- get_coastline(r)
  if (length(v)==0){ plot(r, main=title); if(!is.null(coast)) try(lines(crop(coast, ext(r)), lwd=0.4), silent=TRUE); return() }

  if (palette=="binary"){
    plot(r, main=title, col=c("#f0f0f0","#3b528b"),
         breaks=c(-0.5,0.5,1.5), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else if (palette=="diverging"){
    vmax <- max(abs(range(v, na.rm=TRUE)))
    plot(r, main=title, col=colorRampPalette(c("#3b528b","#f7f7f7","#b40426"))(201),
         zlim=c(-vmax, vmax), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else if (palette=="tri"){
    plot(r, main=title,
         col=c("#b40426","#f7f7f7","#3b528b"), # Loss / Stable / Gain
         breaks=c(0.5,1.5,2.5,3.5), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else {
    plot(r, main=title, col=viridisLite::viridis(200),
         zlim=range(v, na.rm=TRUE), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  }
  if(!is.null(coast)) try(lines(crop(coast, ext(r)), lwd=0.4, col="black"), silent=TRUE)
}

# ------------------------------- Utilities -----------------------------------
write_gt <- function(r, p){ dir.create(dirname(p), TRUE, TRUE); writeRaster(r, p, overwrite=TRUE) }
clamp01  <- function(x){ x <- ifel(x<0,0,x); ifel(x>1,1,x) }
qmask <- function(r, q = 0.75) {
  # prefer global() which works on single-layer rasters
  thr <- tryCatch(
    as.numeric(terra::global(r, fun = quantile, probs = q, na.rm = TRUE)[1, 1]),
    error = function(e) NA_real_
  )
  if (is.na(thr)) {
    # fallback: compute from values() if global() failed for some reason
    v <- terra::values(r, mat = FALSE)
    v <- v[is.finite(v)]
    if (length(v) == 0L) return(r * NA)  # all NA
    thr <- as.numeric(stats::quantile(v, probs = q, na.rm = TRUE, names = FALSE))
  }
  r > thr
}
rmse     <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
mae      <- function(a,b) mean(abs(a-b), na.rm=TRUE)

harmonize <- function(SSF, H2O, SSDM){
  resos <- c(SSF=mean(res(SSF)), H2O=mean(res(H2O)), SSDM=mean(res(SSDM)))
  target <- list(SSF=SSF, H2O=H2O, SSDM=SSDM)[[names(which.min(resos))]]
  conv <- function(r) {
    if (!compareGeom(r, target, stopOnError=FALSE)) {
      r <- project(r, crs(target))
      r <- resample(r, target, method="bilinear")
    }
    r
  }
  stk <- c(conv(SSF), conv(H2O), conv(SSDM)); names(stk) <- c("SSF","H2O","SSDM")
  clamp01(stk)
}

# ------------------------------- Per-run -------------------------------------
analyze_run <- function(ele, run, dirs){
  message(sprintf("Reading %s %s…", ele, run))
  r_ssf  <- rast(p_ssf(ele, run))
  r_h2o  <- rast(p_h2o(ele, run))
  r_ssdm <- rast(p_ssdm(ele, run))

  # Announce inputs (nice for logs)
  message(sprintf("[INPUT] SSF : %s", normalizePath(p_ssf(ele,run),  mustWork=FALSE)))
  message(sprintf("[INPUT] H2O : %s", normalizePath(p_h2o(ele,run),  mustWork=FALSE)))
  message(sprintf("[INPUT] SSDM: %s", normalizePath(p_ssdm(ele,run), mustWork=FALSE)))

  stk <- harmonize(r_ssf, r_h2o, r_ssdm)

  r_mean <- mean(stk, na.rm=TRUE); names(r_mean) <- "mean"
  r_sd   <- stdev(stk, na.rm=TRUE); names(r_sd)   <- "sd"
  r_cv   <- r_sd / (r_mean + 1e-9); names(r_cv)   <- "cv"

  r_alg  <- stk[["H2O"]] - stk[["SSDM"]]; names(r_alg) <- "diff_H2O_minus_SSDM"
  r_sdmM <- (stk[["H2O"]] + stk[["SSDM"]]) / 2
  r_beh  <- r_sdmM - stk[["SSF"]]; names(r_beh) <- "diff_SDMmean_minus_SSF"

  m_ssf  <- qmask(stk[["SSF"]],  opt$qmask)
  m_h2o  <- qmask(stk[["H2O"]],  opt$qmask)
  m_ssdm <- qmask(stk[["SSDM"]], opt$qmask)
  m_ag   <- (m_ssf + m_h2o + m_ssdm) / 3; names(m_ag) <- "agree_prop"

  df <- as.data.frame(stk, xy=FALSE, na.rm=TRUE)

  glb <- tibble(
    elephant = ele, run = run, metric = c("mean","sd","cv"),
    SSF = c(global(stk[["SSF"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["SSF"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["SSF"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["SSF"]],"mean",na.rm=TRUE)[1,1]+1e-9)),
    H2O = c(global(stk[["H2O"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["H2O"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["H2O"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["H2O"]],"mean",na.rm=TRUE)[1,1]+1e-9)),
    SSDM= c(global(stk[["SSDM"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["SSDM"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["SSDM"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["SSDM"]],"mean",na.rm=TRUE)[1,1]+1e-9))
  )

  cor_tbl <- tibble(
    elephant = ele, run = run,
    pair     = c("H2O~SSDM","SSF~H2O","SSF~SSDM","SSF~SDMmean"),
    pearson  = c(cor(df$H2O, df$SSDM, use="pairwise"),
                 cor(df$SSF, df$H2O,  use="pairwise"),
                 cor(df$SSF, df$SSDM, use="pairwise"),
                 cor(df$SSF, (df$H2O+df$SSDM)/2, use="pairwise")),
    spearman = c(cor(df$H2O, df$SSDM, use="pairwise", method="spearman"),
                 cor(df$SSF, df$H2O,  use="pairwise", method="spearman"),
                 cor(df$SSF, df$SSDM, use="pairwise", method="spearman"),
                 cor(df$SSF, (df$H2O+df$SSDM)/2, use="pairwise", method="spearman")),
    rmse     = c(rmse(df$H2O, df$SSDM),
                 rmse(df$SSF, df$H2O),
                 rmse(df$SSF, df$SSDM),
                 rmse(df$SSF, (df$H2O+df$SSDM)/2)),
    mae      = c(mae(df$H2O, df$SSDM),
                 mae(df$SSF, df$H2O),
                 mae(df$SSF, df$SSDM),
                 mae(df$SSF, (df$H2O+df$SSDM)/2))
  )

  # Pairwise Jaccard@Q
  jac <- {
    inter <- global((m_ssf & m_h2o), "sum", na.rm=TRUE)[1,1]; uni <- global((m_ssf | m_h2o), "sum", na.rm=TRUE)[1,1]
    j1 <- ifelse(uni==0, NA_real_, inter/uni)
    inter <- global((m_ssf & m_ssdm), "sum", na.rm=TRUE)[1,1]; uni <- global((m_ssf | m_ssdm), "sum", na.rm=TRUE)[1,1]
    j2 <- ifelse(uni==0, NA_real_, inter/uni)
    inter <- global((m_h2o & m_ssdm), "sum", na.rm=TRUE)[1,1]; uni <- global((m_h2o | m_ssdm), "sum", na.rm=TRUE)[1,1]
    j3 <- ifelse(uni==0, NA_real_, inter/uni)
    tibble(elephant=ele, run=run, q=opt$qmask,
           pair=c("SSF~H2O","SSF~SSDM","H2O~SSDM"),
           jaccard_q=c(j1,j2,j3))
  }

  # Write rasters
  write_gt(r_mean, file.path(dirs$combined, sprintf("%s_%s_mean.tif", ele, run)))
  write_gt(r_sd,   file.path(dirs$combined, sprintf("%s_%s_sd.tif",   ele, run)))
  write_gt(r_cv,   file.path(dirs$combined, sprintf("%s_%s_cv.tif",   ele, run)))
  write_gt(r_alg,  file.path(dirs$combined, sprintf("%s_%s_diff_H2O_minus_SSDM.tif", ele, run)))
  write_gt(r_beh,  file.path(dirs$combined, sprintf("%s_%s_diff_SDMmean_minus_SSF.tif", ele, run)))
  write_gt(m_ag,   file.path(dirs$combined, sprintf("%s_%s_agreement_prop_q%02d.tif", ele, run, round(opt$qmask*100))))

  # PNGs
  if (!opt$skip_figs){
    save_raster_map(r_mean, sprintf("%s %s — Mean", ele, run), file.path(dirs$combined, sprintf("%s_%s_mean.png", ele, run)))
    save_raster_map(r_sd,   sprintf("%s %s — SD",   ele, run), file.path(dirs$combined, sprintf("%s_%s_sd.png",   ele, run)))
    save_raster_map(r_cv,   sprintf("%s %s — CV",   ele, run), file.path(dirs$combined, sprintf("%s_%s_cv.png",   ele, run)))
    save_raster_map(m_ag,   sprintf("%s %s — Agreement (Q%02d)", ele, run, round(opt$qmask*100)),
                    file.path(dirs$combined, sprintf("%s_%s_agreement_q%02d.png", ele, run, round(opt$qmask*100))))
    save_raster_map(r_alg,  sprintf("%s %s — H2O − SSDM", ele, run),
                    file.path(dirs$combined, sprintf("%s_%s_diff_H2O_minus_SSDM.png", ele, run)), "diverging")
    save_raster_map(r_beh,  sprintf("%s %s — SDM̄ − SSF", ele, run),
                    file.path(dirs$combined, sprintf("%s_%s_diff_SDMmean_minus_SSF.png", ele, run)), "diverging")
  }

  list(stack=stk, summaries=glb, cors=cor_tbl, jacc=jac)
}

# ----------------------------- Temporal A vs B --------------------------------
temporal_compare <- function(ele, Astk, Bstk, dirs){
  methods <- c("SSF","H2O","SSDM")
  jac <- list()
  for (m in methods){
    rA <- Astk[[m]]; rB <- Bstk[[m]]
    d  <- rB - rA
    write_gt(d, file.path(dirs[[m]], sprintf("%s_delta_BminusA.tif", ele)))
    if (!opt$skip_figs)
      save_raster_map(d, sprintf("%s Δ (B−A) — %s", ele, m),
                      file.path(dirs[[m]], sprintf("%s_delta_BminusA.png", ele)), "diverging")

    mA <- qmask(rA, opt$qmask); mB <- qmask(rB, opt$qmask)
    inter <- global(mA & mB, "sum", na.rm=TRUE)[1,1]
    uni   <- global(mA | mB, "sum", na.rm=TRUE)[1,1]
    jac[[m]] <- tibble(elephant=ele, method=m, q=opt$qmask, jaccard_AvB=ifelse(uni==0, NA_real_, inter/uni))

    # Change classes: 1=Loss, 2=Stable, 3=Gain
    loss   <- (mA==1) & (mB==0)
    stable <- (mA==1) & (mB==1)
    gain   <- (mA==0) & (mB==1)
    change <- (loss*1) + (stable*2) + (gain*3)
    names(change) <- paste0("change_", m, "_Q", round(opt$qmask*100))
    write_gt(change, file.path(dirs$combined, sprintf("%s_change_%s_q%02d.tif", ele, m, round(opt$qmask*100))))
    if (!opt$skip_figs)
      save_raster_map(change, sprintf("%s Change classes — %s (Q%02d)", ele, m, round(opt$qmask*100)),
                      file.path(dirs$combined, sprintf("%s_change_%s_q%02d.png", ele, m, round(opt$qmask*100))), "tri")
  }
  bind_rows(jac)
}

# -------------------------- DEBUG: show resolved paths ------------------------
debug_inputs <- function(ele, run) {
  h <- p_h2o(ele, run); s <- p_ssdm(ele, run); f <- p_ssf(ele, run)
  cat("\n[DEBUG]", ele, run,
      "\n H2O :",  normalizePath(h, mustWork = FALSE), "exists?", file.exists(h),
      "\n SSDM:",  normalizePath(s, mustWork = FALSE), "exists?", file.exists(s),
      "\n SSF :",  normalizePath(f, mustWork = FALSE), "exists?", file.exists(f), "\n")
}

# --------------------------------- MAIN ---------------------------------------
# Parse elephant set and normalize to base codes (E3, E4, ...)
elephants <- toupper(gsub("[^A-Z0-9,]", "", opt$elephant))
if (identical(elephants, "ALL")) {
  elephants <- c("E1","E2","E3","E4","E5","E6")
} else {
  elephants <- strsplit(elephants, ",", fixed = TRUE)[[1]]
}
elephants <- elephants[nzchar(elephants)]

# If user passed E3A/E4B, strip the trailing run letter (avoid E3AA/E4BB)
normalize_ele <- function(s) {
  s <- trimws(s)
  if (nchar(s) >= 2) {
    last <- substr(s, nchar(s), nchar(s))
    if (last %in% c("A","B")) return(substr(s, 1, nchar(s)-1))
  }
  s
}
elephants <- unique(vapply(elephants, normalize_ele, FUN.VALUE = character(1)))

for (ele in elephants) {
  dirs <- mk_outdirs(ele)

  message(sprintf("\n== Uncertainty decomposition: %s ==", ele))

  # explicit debug: show exactly what paths will be used
  debug_inputs(ele, runA)
  debug_inputs(ele, runB)

  hasA <- exists_all(ele, runA)
  hasB <- exists_all(ele, runB)

  if (!hasA && !hasB) {
    warning(sprintf("Skipping %s: neither %s nor %s inputs found.", ele, runA, runB))
    next
  }

  # ----- Per-run analysis -----
  A <- NULL; B <- NULL
  if (hasA) {
    A <- analyze_run(ele, runA, dirs)
  } else {
    message(sprintf("No %s run for %s — skipping A.", runA, ele))
  }

  if (hasB) {
    B <- analyze_run(ele, runB, dirs)
  } else {
    message(sprintf("No %s run for %s — skipping B.", runB, ele))
  }

  # ----- Tables per run -----
  tbl_dir <- file.path(opt$outdir, ele, "tables"); dir.create(tbl_dir, TRUE, TRUE)
  if (!is.null(A)) {
    write_csv(A$summaries, file.path(tbl_dir, sprintf("%s_global_stats_%s.csv", ele, runA)))
    write_csv(A$cors,      file.path(tbl_dir, sprintf("%s_pairwise_correlations_%s.csv", ele, runA)))
    write_csv(A$jacc,      file.path(tbl_dir, sprintf("%s_jaccard_q%02d_%s.csv", ele, round(opt$qmask*100), runA)))
  }
  if (!is.null(B)) {
    write_csv(B$summaries, file.path(tbl_dir, sprintf("%s_global_stats_%s.csv", ele, runB)))
    write_csv(B$cors,      file.path(tbl_dir, sprintf("%s_pairwise_correlations_%s.csv", ele, runB)))
    write_csv(B$jacc,      file.path(tbl_dir, sprintf("%s_jaccard_q%02d_%s.csv", ele, round(opt$qmask*100), runB)))
  }

  # ----- Temporal only if both runs exist -----
  if (!is.null(A) && !is.null(B)) {
    TEMP <- temporal_compare(ele, A$stack, B$stack, dirs)
    write_csv(TEMP, file.path(tbl_dir, sprintf("%s_temporal_jaccard_q%02d.csv", ele, round(opt$qmask*100))))
  } else {
    message(sprintf("Temporal skipped for %s (need both %s and %s).", ele, runA, runB))
  }

  message(sprintf("Finished %s.", ele))
}

message("\nAll done.\n")
