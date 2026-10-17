# Shared validation and summaries for DATA309_EDA.R; no files are changed here.
require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

validate_cells <- function(data, label) {
  require_columns(data, c("publicid", "cell_id", "origintime", "event_longitude",
    "event_latitude", "event_depth_km", "report_longitude", "report_latitude",
    "reported_mmi", "report_count", "Mw", "Repi_km", "Rhypo_km"), label)
  if (!nrow(data)) stop(label, " has no observations")
  if (anyNA(data[c("publicid", "cell_id")]) ||
      anyDuplicated(data[c("publicid", "cell_id")])) {
    stop(label, " must have unique, non-missing earthquake-cell keys")
  }
  numeric_fields <- c("event_longitude", "event_latitude", "event_depth_km",
    "report_longitude", "report_latitude", "reported_mmi", "report_count", "Mw",
    "Repi_km", "Rhypo_km", intersect("Vs30_m_s", names(data)))
  for (field in numeric_fields) {
    x <- data[[field]]
    if (!is.numeric(x) || any(!is.finite(x) & !is.na(x))) stop("Invalid ", field)
  }
  if (!inherits(data$origintime, "POSIXt") || anyNA(data$origintime)) {
    stop(label, ": origintime must contain valid timestamps")
  }
  if (any(data$report_count < 1 | data$report_count != floor(data$report_count), na.rm = TRUE)) {
    stop(label, ": report_count must contain positive integers")
  }
  # Aggregated medians can be fractional, even though the current file has integers.
  if (any(data$reported_mmi <= 0 | data$reported_mmi > 12, na.rm = TRUE)) stop("Invalid MMI")
  for (field in c("event_depth_km", "Repi_km", "Rhypo_km")) {
    if (any(data[[field]] < 0, na.rm = TRUE)) stop("Negative ", field)
  }
  for (field in c("event_longitude", "report_longitude")) {
    if (any(abs(data[[field]]) > 180, na.rm = TRUE)) stop("Invalid longitude")
  }
  for (field in c("event_latitude", "report_latitude")) {
    if (any(abs(data[[field]]) > 90, na.rm = TRUE)) stop("Invalid latitude")
  }
  if ("Vs30_m_s" %in% names(data) && any(data$Vs30_m_s <= 0, na.rm = TRUE)) stop("Invalid Vs30")
  event_fields <- c("origintime", "event_longitude", "event_latitude", "event_depth_km", "Mw")
  event_values <- unique(data[c("publicid", event_fields)])
  if (anyDuplicated(event_values$publicid)) stop(label, ": inconsistent earthquake attributes")
  invisible(TRUE)
}

validate_split_ids <- function(train, test) {
  for (ids in list(train, test)) {
    if (!is.character(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) {
      stop("Saved earthquake IDs must be unique, non-missing character vectors")
    }
  }
  if (length(intersect(train, test))) stop("Training and test earthquakes overlap")
  invisible(TRUE)
}

numeric_summary <- function(x) {
  observed <- x[is.finite(x)]
  if (!length(observed)) return(data.frame(n = 0, missing = length(x), min = NA_real_,
    q25 = NA_real_, median = NA_real_, mean = NA_real_, q75 = NA_real_, max = NA_real_))
  q <- quantile(observed, c(0, .25, .5, .75, 1), names = FALSE)
  data.frame(n = length(observed), missing = sum(!is.finite(x)), min = q[1],
    q25 = q[2], median = q[3], mean = mean(observed), q75 = q[4], max = q[5])
}

safe_spearman <- function(x, y, minimum_n = 10L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < minimum_n || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = "spearman")
}

# Publish only after all calculations succeed. Refuse an unknown directory so
# unrelated user files cannot be replaced. A failed rename restores the old run.
publish_eda <- function(stage, destination) {
  marker <- ".DATA309_EDA"
  if (dir.exists(destination)) {
    if (!file.exists(file.path(destination, marker))) {
      stop("Output directory is not owned by DATA309_EDA; choose a new --output-dir")
    }
    owned <- read.csv(file.path(destination, "output_manifest.csv"), stringsAsFactors = FALSE)$file
    extra <- setdiff(list.files(destination, all.files = TRUE, no.. = TRUE),
      c(owned, marker, "output_manifest.csv"))
    if (length(extra)) stop("Output directory contains untracked files: ", paste(extra, collapse = ", "))
  }
  backup <- tempfile("eda-backup-", tmpdir = dirname(destination))
  existed <- dir.exists(destination)
  if (existed && !file.rename(destination, backup)) stop("Could not back up previous EDA")
  if (!file.rename(stage, destination)) {
    if (existed) file.rename(backup, destination)
    stop("Could not publish EDA; previous output restored")
  }
  if (existed) unlink(backup, recursive = TRUE)
  invisible(destination)
}
