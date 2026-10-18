# 01_combine_geonet_catalogue
# Thomas Rautao 2026

# Combine the monthly GeoNet files into one earthquake catalogue

library(tidyverse)
library(arrow)

input_dir <- file.path("data_raw", "geonet_catalogue")
output_file <- file.path("data_raw", "geonet_catalogue.parquet")

# Sorting the filenames sorts the monthly files chronologically because 
# they are formatted YYYY_MM.csv.
files <- list.files(
  input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
) |>
  sort()

if (length(files) == 0) {
  stop("No monthly GeoNet catalogue files were found. Run Script 00 first.")
}

# columns are initially read as text to prevent row binding from failing.
# Required fields are converted after all files are joined.
read_catalogue_file <- function(path) {
  read_csv(
    path,
    col_types = cols(.default = col_character()),
    na = c("", "NA", "N/A", "NaN", "NULL", "null"),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    rename_with(~ str_to_lower(str_trim(.x)))
}

catalogue_raw <- map_dfr(files, read_catalogue_file)

required_columns <- c(
  "publicid", "origintime", "latitude", "longitude", "depth"
)
missing_columns <- setdiff(required_columns, names(catalogue_raw))

if (length(missing_columns) > 0) {
  stop("Missing catalogue columns: ", paste(missing_columns, collapse = ", "))
}

catalogue <- catalogue_raw |>
  transmute(
    publicid = na_if(str_trim(publicid), ""),
    origintime = lubridate::ymd_hms(origintime, tz = "UTC", quiet = TRUE),
    latitude = suppressWarnings(parse_double(latitude)),
    longitude = suppressWarnings(parse_double(longitude)),
    depth = suppressWarnings(parse_double(depth))
  ) |>
  filter(!is.na(publicid)) |>
  distinct(publicid, .keep_all = TRUE) |>
  arrange(origintime)

if (nrow(catalogue) == 0) {
  stop("The monthly files did not contain any usable earthquakes.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(catalogue, output_file)

message("Saved ", nrow(catalogue), " earthquakes to ", output_file)
