# 04_build_nz_fr_iar_dataset
# Thomas Rautao 2026

# Build the earthquake-by-Felt-RAPID-cell dataset used for analysis.

# Each row links one earthquake to one Felt RAPID geographic cell. 

library(tidyverse)
library(arrow)

analysis_start_date <- as.POSIXct("2016-09-01 00:00:00", tz = "UTC")
minimum_moment_magnitude <- 3.5

event_file <- file.path(
  "data_processed",
  "geonet_catalogue_with_mw.parquet"
)
felt_file <- file.path("data_raw", "felt_rapid_observations.parquet")
output_file <- file.path("data_processed", "NZ_FR_IAR_dataset.parquet")

if (!file.exists(event_file)) {
  stop("Moment-tensor catalogue not found. Run Script 03 first.")
}

if (!file.exists(felt_file)) {
  stop("Felt RAPID observations not found. Run Script 02 first.")
}

events <- read_parquet(event_file) |>
  as_tibble()

felt <- read_parquet(felt_file) |>
  as_tibble()

required_event_columns <- c(
  "publicid", "origintime", "latitude", "longitude", "depth", "Mw"
)
required_felt_columns <- c(
  "publicid", "cell_id", "report_longitude", "report_latitude",
  "reported_mmi", "report_count"
)

missing_event_columns <- setdiff(required_event_columns, names(events))
missing_felt_columns <- setdiff(required_felt_columns, names(felt))

if (length(missing_event_columns) > 0) {
  stop("Missing event columns: ", paste(missing_event_columns, collapse = ", "))
}

if (length(missing_felt_columns) > 0) {
  stop("Missing Felt RAPID columns: ", paste(missing_felt_columns, collapse = ", "))
}

# Rename earthquake coordinates to distinguish them from the coordinates
# representing the centre of each Felt RAPID reporting cell.
events <- events |>
  transmute(
    publicid,
    origintime,
    event_longitude = suppressWarnings(as.numeric(longitude)),
    event_latitude = suppressWarnings(as.numeric(latitude)),
    event_depth_km = suppressWarnings(as.numeric(depth)),
    Mw = suppressWarnings(as.numeric(Mw))
  )

felt <- felt |>
  transmute(
    publicid,
    cell_id,
    report_longitude = suppressWarnings(as.numeric(report_longitude)),
    report_latitude = suppressWarnings(as.numeric(report_latitude)),
    reported_mmi = suppressWarnings(as.numeric(reported_mmi)),
    report_count = suppressWarnings(as.integer(report_count))
  ) |>
  distinct(publicid, cell_id, .keep_all = TRUE)

# remove missing event matches rows
iar <- felt |>
  left_join(events, by = "publicid") |>
  mutate(
    valid_report_coordinates =
      !is.na(report_longitude) &
      !is.na(report_latitude) &
      between(report_longitude, -180, 180) &
      between(report_latitude, -90, 90),
    valid_event_coordinates =
      !is.na(event_longitude) &
      !is.na(event_latitude) &
      between(event_longitude, -180, 180) &
      between(event_latitude, -90, 90),
    valid_event_depth =
      !is.na(event_depth_km) &
      between(event_depth_km, 0, 700)
  )

# S2 geometry accounts for the curvature of Earth and handles the NZ 
# longitude range 
iar$Repi_km <- NA_real_
distance_rows <- which(
  iar$valid_report_coordinates & iar$valid_event_coordinates
)

if (length(distance_rows) > 0) {
  sf::sf_use_s2(TRUE)

  report_points <- sf::st_as_sf(
    iar[distance_rows, ],
    coords = c("report_longitude", "report_latitude"),
    crs = 4326,
    remove = FALSE
  )

  event_points <- sf::st_as_sf(
    iar[distance_rows, ],
    coords = c("event_longitude", "event_latitude"),
    crs = 4326,
    remove = FALSE
  )

  iar$Repi_km[distance_rows] <- sf::st_distance(
    report_points,
    event_points,
    by_element = TRUE
  ) |>
    units::set_units("km") |>
    units::drop_units() |>
    as.numeric()
}

# Straight-line distance from earthquake hypocentre to the Felt cell
# Rhypo = sqrt(Repi^2 + depth^2)
iar <- iar |>
  mutate(
    Rhypo_km = if_else(
      !is.na(Repi_km) & valid_event_depth,
      sqrt(Repi_km^2 + event_depth_km^2),
      NA_real_
    )
  ) |>
  filter(
    !is.na(origintime),
    origintime >= analysis_start_date,
    !is.na(Mw),
    Mw > minimum_moment_magnitude,
    !is.na(reported_mmi),
    between(reported_mmi, 1, 12),
    !is.na(report_count),
    report_count > 0,
    valid_report_coordinates,
    valid_event_coordinates,
    valid_event_depth,
    !is.na(Repi_km),
    !is.na(Rhypo_km)
  ) |>
  select(
    publicid,
    cell_id,
    origintime,
    event_longitude,
    event_latitude,
    event_depth_km,
    report_longitude,
    report_latitude,
    reported_mmi,
    report_count,
    Mw,
    Repi_km,
    Rhypo_km
  ) |>
  arrange(publicid, Repi_km, cell_id)

if (nrow(iar) == 0) {
  stop("The analysis dataset contains no eligible observations.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(iar, output_file)

message(
  "Saved ", nrow(iar), " Felt RAPID cells from ",
  n_distinct(iar$publicid), " earthquakes"
)
