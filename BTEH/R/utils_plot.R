suppressPackageStartupMessages({
  library(ggplot2); library(scales)
})

theme_pub <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b=6)),
      plot.subtitle = element_text(margin = margin(b=8)),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(1.0, "cm"),    # slightly wider
      legend.key.height = unit(0.6, "cm")    # slightly taller
    )
}

# Simple helper: annotate values on bars
geom_bar_labels <- function(fmt = label_number(accuracy = 0.01), vjust = -0.2, size = 3) {
  geom_text(aes(label = fmt(after_stat(y))), stat = "sum", vjust = vjust, size = size)
}
