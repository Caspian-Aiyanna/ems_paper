suppressPackageStartupMessages({
  library(fs)
  library(readr)
})

dir_ensure <- function(path) {
  if (!fs::dir_exists(path)) fs::dir_create(path, recurse = TRUE)
  invisible(path)
}

path_proj <- function(...) {
  p <- fs::path(...)
  fs::path_norm(p)
}

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_write_csv <- function(df, path) {
  dir_ensure(fs::path_dir(path))
  readr::write_csv(df, path)
  invisible(path)
}

log_line <- function(msg, log_file = NULL) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg)
  message(line)
  if (!is.null(log_file)) {
    dir_ensure(fs::path_dir(log_file))
    cat(line, file = log_file, sep = "\n", append = TRUE)
  }
  invisible(line)
}
