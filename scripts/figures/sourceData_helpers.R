############################################################
# Script: sourceData_helpers.R
# Purpose: Shared helper for exporting the tabular data behind
#          each figure panel/table to source_data/, for the
#          journal-required Source Data workbook. Sourced by
#          every script in scripts/figures/.
############################################################

#' Write a tidy data frame as a source-data CSV for one figure panel/table.
#'
#' @param df      data.frame/tibble to export (caller pre-selects/renames cols)
#' @param name    short id, e.g. "AMMI_fig4_F" -> source_data/main_figs/AMMI_fig4_F.csv
#' @param section "main_figs" or "supp_figs"
save_source_data <- function(df, name, section = c("main_figs", "supp_figs")) {
  section <- match.arg(section)
  dir <- file.path("source_data", section)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  df <- as.data.frame(df)
  # defensive: flatten any list-columns (e.g. leftover nested/boot list-cols)
  list_cols <- vapply(df, is.list, logical(1))
  if (any(list_cols)) {
    df[list_cols] <- lapply(df[list_cols], function(x) vapply(x, toString, character(1)))
  }

  out <- file.path(dir, paste0(name, ".csv"))
  write.csv(df, out, row.names = FALSE)
  invisible(out)
}

#' Five-number summary matching what geom_boxplot() actually draws.
#'
#' Hinges via quantile() (type 7, ggplot2's default); whiskers extended to
#' the most extreme non-outlier value within 1.5*IQR of the hinges
#' (ggplot2's default coef). NOT simple min/max when outliers are present.
#'
#' @param y numeric vector of the values the boxplot is drawn on (i.e. the
#'   already-transformed plotted values, e.g. -log10(p), not the raw values,
#'   if the plot itself transforms them)
box_stats <- function(y) {
  qs  <- stats::quantile(y, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  iqr <- qs[4] - qs[2]
  lower_fence <- qs[2] - 1.5 * iqr
  upper_fence <- qs[4] + 1.5 * iqr
  in_range <- y >= lower_fence & y <= upper_fence
  data.frame(
    n             = length(y),
    min           = qs[1],
    q1            = qs[2],
    median        = qs[3],
    q3            = qs[4],
    max           = qs[5],
    whisker_lower = if (any(in_range)) min(y[in_range], na.rm = TRUE) else qs[2],
    whisker_upper = if (any(in_range)) max(y[in_range], na.rm = TRUE) else qs[4],
    n_outliers    = sum(!in_range)
  )
}
