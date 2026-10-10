# 00_get_geonet_catalogue
# Thomas Rautao 2026

# Download monthly GeoNet earthquake catalogue files

library(readr)
library(lubridate)

output_dir <- file.path("data_raw", "geonet_catalogue")

# Request filter
minimum_catalogue_magnitude <- 3.0

# NZ bounding box, west, south, east, north.
# Coordinates taken from the GeoNet QuakeSearch page
bbox <- "163.5205,-49.1817,-176.9238,-32.2871"

# Default False, set True if need to redo all months after changing parameters
overwrite_existing <- FALSE

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Felt RAPID reports from 1 September 2016, as per project outline
current_month <- floor_date(Sys.Date(), "month")
months <- seq(
  as.Date("2016-09-01"),
  current_month,
  by = "month"
)

for (i in seq_along(months)) {
  start_date <- months[i]
  end_date <- if (i < length(months)) months[i + 1] else Sys.Date() + 1
  output_file <- file.path(
    output_dir,
    paste0(format(start_date, "%Y_%m"), ".csv")
  )

  if (
    file.exists(output_file) &&
      !overwrite_existing &&
      start_date < current_month
  ) {
    message("Using existing file: ", output_file)
    next
  }

  message("Downloading GeoNet catalogue: ", start_date, " to ", end_date)

  url <- paste0(
    "https://quakesearch.geonet.org.nz/csv?",
    "bbox=", bbox,
    "&minmag=", minimum_catalogue_magnitude,
    "&startdate=", start_date, "T00:00:00",
    "&enddate=", end_date, "T00:00:00"
  )

  catalogue_month <- read_csv(url, show_col_types = FALSE)
  write_csv(catalogue_month, output_file)

  message("Saved ", nrow(catalogue_month), " earthquakes to ", output_file)

  Sys.sleep(1)
}
