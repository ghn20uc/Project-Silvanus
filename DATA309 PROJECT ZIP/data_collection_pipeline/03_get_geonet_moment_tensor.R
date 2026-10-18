# 03_get_geonet_moment_tensor
# Thomas Rautao 2026

# Add GeoNet moment-tensor Mw values to the earthquake catalogue

# Mw is the only used magnitude (not the earthquake catalogue magnitude)

library(tidyverse)
library(arrow)

input_file <- file.path("data_raw", "geonet_catalogue.parquet")
output_file <- file.path(
  "data_processed",
  "geonet_catalogue_with_mw.parquet"
)

# Keep the original source, including solutions that are later excluded or
# superseded. Script 06 needs this to count source rows, rather than counting
# only earthquakes that matched the study catalogue.
mt_source_file <- file.path(
  "data_raw", "moment_tensor", "GeoNet_CMT_solutions.csv"
)
mt_metadata_file <- file.path(
  "data_raw", "moment_tensor", "source_metadata.csv"
)

mt_url <- paste0(
  "https://raw.githubusercontent.com/GeoNet/data/main/moment-tensor/GeoNet_CMT_solutions.csv"
)

if (!file.exists(input_file)) {
  stop("GeoNet catalogue not found. Run Script 01 first.")
}

catalogue <- read_parquet(input_file) |>
  as_tibble()

message("Downloading GeoNet moment-tensor catalogue")

# Download to a temporary file first. A failed download shouldn't overwrite
# an already saved source. Read and validate download before saving it.
dir.create(dirname(mt_source_file), recursive = TRUE, showWarnings = FALSE)
mt_download_file <- tempfile(
  pattern = "GeoNet_CMT_", tmpdir = dirname(mt_source_file), fileext = ".csv"
)
download.file(mt_url, mt_download_file, mode = "wb", quiet = TRUE)
mt_downloaded_at_utc <- format(
  Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
)

# All source columns are read as text so missing-value and numeric
# conversion can be done below.
mt_raw <- read_csv(
  mt_download_file,
  col_types = cols(.default = col_character()),
  na = c("", "NA", "N/A", "NaN", "NULL", "null", "-"),
  show_col_types = FALSE,
  progress = FALSE
) |>
  rename_with(~ str_to_lower(str_trim(.x)))

# Standardising naming from event_id/publicid naming
if ("event_id" %in% names(mt_raw) && !"publicid" %in% names(mt_raw)) {
  mt_raw <- rename(mt_raw, publicid = event_id)
}

# Mw is the analysis magnitude. Variance reduction and station count are used
# to choose between multiple valid moment-tensor solutions for the same
# earthquake 
required_columns <- c("publicid", "mw", "vr", "ns")
missing_columns <- setdiff(required_columns, names(mt_raw))

if (length(missing_columns) > 0) {
  stop("Missing moment-tensor columns: ", paste(missing_columns, collapse = ", "))
}

mt_candidates <- mt_raw |>
  transmute(
    source_row = row_number(),
    publicid = na_if(str_trim(publicid), ""),
    Mw = suppressWarnings(parse_double(mw)),
    variance_reduction = suppressWarnings(parse_double(vr)),
    station_count = suppressWarnings(parse_double(ns))
  ) |>
  filter(
    !is.na(publicid),
    publicid != "9999999",
    !is.na(Mw),
    between(Mw, 0, 10)
  )

if (nrow(mt_candidates) == 0) {
  stop("No usable moment-tensor values were found.")
}

# When GeoNet supplies more than one CMT solution for an earthquake, the
# solution with the highest variance reduction is selected because it provides
# the best fit. Station count is used to break a tie. source row provides 
# a second tie-break if needed.
moment_magnitudes <- mt_candidates |>
  arrange(
    publicid,
    desc(coalesce(variance_reduction, -Inf)),
    desc(coalesce(station_count, -Inf)),
    desc(source_row)
  ) |>
  distinct(publicid, .keep_all = TRUE) |>
  select(publicid, Mw) |>
  arrange(publicid)

# Retain every downloaded catalogue earthquake. Events without a
# moment-tensor solution receive Mw = NA.
catalogue_with_mw <- catalogue |>
  left_join(moment_magnitudes, by = "publicid")

matched_count <- sum(!is.na(catalogue_with_mw$Mw))

if (matched_count == 0) {
  stop("No catalogue earthquakes matched a moment-tensor value.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(catalogue_with_mw, output_file)

# Save the exact downloaded bytes, rather than writing the filtered or
# renamed table back to CSV,preserves every original source column.
if (!file.copy(mt_download_file, mt_source_file, overwrite = TRUE)) {
  stop("Could not save the original moment-tensor CSV.")
}
unlink(mt_download_file)

# Save diagnostic data
write_csv(
  tibble(
    source_url = mt_url,
    downloaded_at_utc = mt_downloaded_at_utc,
    raw_file = mt_source_file,
    source_rows = nrow(mt_raw),
    selected_event_count = nrow(moment_magnitudes),
    matched_event_count = matched_count,
    raw_file_md5 = unname(tools::md5sum(mt_source_file)),
    catalogue_with_mw_file_md5 = unname(tools::md5sum(output_file))
  ),
  mt_metadata_file,
  na = "NA"
)

message(
  "Matched moment-tensor Mw to ",
  matched_count,
  " of ", nrow(catalogue_with_mw),
  " catalogue earthquakes"
)
