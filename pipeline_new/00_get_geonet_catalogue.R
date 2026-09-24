# Download monthly GeoNet earthquake catalogue files.

library(readr)
library(lubridate)

output_dir <- file.path("data_raw", "geonet_catalogue")
minimum_catalogue_magnitude <- 3.5
bbox <- "163.5205,-49.1817,-176.9238,-32.2871"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

months <- seq(
  as.Date("2016-09-01"),
  floor_date(Sys.Date(), "month"),
  by = "month"
)

for (i in seq_along(months)) {
  start_date <- months[i]
  end_date <- if (i < length(months)) months[i + 1] else Sys.Date() + 1
  output_file <- file.path(output_dir, paste0(format(start_date, "%Y_%m"), ".csv"))

  if (file.exists(output_file)) {
    message("Skipping ", output_file)
    next
  }

  message("Downloading ", start_date, " to ", end_date)

  url <- paste0(
    "https://quakesearch.geonet.org.nz/csv?",
    "bbox=", bbox,
    "&minmag=", minimum_catalogue_magnitude,
    "&startdate=", start_date, "T00:00:00",
    "&enddate=", end_date, "T00:00:00"
  )

  catalogue_month <- read_csv(url, show_col_types = FALSE)
  write_csv(catalogue_month, output_file)

  message("Downloaded ", nrow(catalogue_month), " earthquakes")
  Sys.sleep(1)
}
