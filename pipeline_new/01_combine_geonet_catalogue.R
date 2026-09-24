# Combine the monthly GeoNet catalogue files.

library(tidyverse)
library(arrow)

input_dir <- file.path("data_raw", "geonet_catalogue")
output_file <- file.path("data_raw", "geonet_catalogue.parquet")

files <- list.files(
  input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
) |>
  sort()

if (length(files) == 0) {
  stop("No monthly GeoNet catalogue files were found in ", input_dir, ".")
}

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
  "publicid", "origintime", "latitude", "longitude", "depth", "magnitude"
)
missing_columns <- setdiff(required_columns, names(catalogue_raw))

if (length(missing_columns) > 0) {
  stop("Missing catalogue columns: ", paste(missing_columns, collapse = ", "))
}

catalogue <- catalogue_raw |>
  mutate(
    publicid = na_if(str_trim(publicid), ""),
    origintime = ymd_hms(origintime, tz = "UTC", quiet = TRUE),
    across(
      c(latitude, longitude, depth, magnitude),
      ~ suppressWarnings(parse_double(.x))
    )
  ) |>
  filter(!is.na(publicid)) |>
  distinct(publicid, .keep_all = TRUE) |>
  mutate(year = year(origintime))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(catalogue, output_file)

message("Saved ", format(nrow(catalogue), big.mark = ","), " earthquakes to ", output_file)
