# Download monthly GeoNet earthquake catalogue files

# NOTES TO SELF:
# minimum-magnitude is set to 5.0, lowering this exponentially increases the time for the script to run
# the bounding box is set to cover all of NZ (South, North, Stewart - "163.5205,-49.1817,-176.9238,-32.2871"), 
# but crosses a date-line boundary, so worth checking if issues arise

library(readr)
library(lubridate)


# Download settings

output_dir <- file.path(
  "data_raw",
  "geonet_catalogue"
)

minimum_magnitude <- 3.0

# Set the bounding box region
bbox <- "163.5205,-49.1817,-176.9238,-32.2871"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Monthly dates

months <- seq(
  as.Date("2016-11-01"),
  floor_date(Sys.Date(), "month"),
  by = "month"
)


# Download files for every month

for (i in seq_along(months)) {
  
  start <- months[i]
  
  if (i < length(months)) {
    end <- months[i + 1]
  } else {
    end <- Sys.Date() + 1
  }
  
  
  outfile <- file.path(
    output_dir,
    paste0(
      format(start, "%Y_%m"),
      ".csv"
    )
  )
  
  
  if (file.exists(outfile)) {
    message("Skipping ", outfile)
    next
  }
  
  
  message(
    "Downloading ",
    start,
    " to ",
    end
  )
  
  
  url <- paste0(
    "https://quakesearch.geonet.org.nz/csv?",
    "bbox=",
    bbox,
    "&minmag=",
    minimum_magnitude,
    "&startdate=",
    start,
    "T00:00:00",
    "&enddate=",
    end,
    "T00:00:00"
  )
  
  
  dat <- read_csv(
    url,
    show_col_types = FALSE
  )
  
  
  write_csv(
    dat,
    outfile
  )
  
  
  message(
    "Downloaded ",
    nrow(dat),
    " earthquakes"
  )
  
  
  Sys.sleep(1)
  
}