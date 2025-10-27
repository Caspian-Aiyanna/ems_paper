suppressPackageStartupMessages({
  library(yaml)
})

read_config <- function(path = "config.yml") {
  yaml::read_yaml(path)
}

set_mode <- function(cfg) {
  mode <- toupper(cfg$mode %||% "REPRO")
  if (mode == "REPRO") {
    # keep it single-threaded / deterministic where applicable
    Sys.setenv("OMP_NUM_THREADS" = "1",
               "MKL_NUM_THREADS" = "1",
               "OPENBLAS_NUM_THREADS" = "1")
  }
  invisible(mode)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

seed_for <- function(tag, base = 1L) {
  # deterministic mapping of a tag to a seed
  h <- as.integer(abs(stats::setNames(utils::adist(tag, tag), NULL))) + as.integer(base)
  set.seed(h)
  invisible(h)
}

save_or_load <- function(path, value_expr) {
  if (file.exists(path)) {
    return(readr::read_rds(path))
  } else {
    val <- eval(substitute(value_expr))
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    readr::write_rds(val, path)
    return(val)
  }
}

skip_if_done <- function(dir, files = character()) {
  all(file.exists(file.path(dir, files)))
}
