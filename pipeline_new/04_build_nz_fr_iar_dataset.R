# Build the earthquake-by-Felt-RAPID-cell dataset used for modelling.

library(tidyverse)
library(arrow)

analysis_start_date <- as.POSIXct("2016-09-01 00:00:00", tz = "UTC")
analysis_end_date <- as.POSIXct(NA_character_, tz = "UTC")
minimum_moment_magnitude <- 3.5
maximum_epicentral_distance_km <- Inf
minimum_reports_per_cell <- 1L

event_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt.parquet"
)
felt_file <- file.path("data_raw", "felt_rapid_observations.parquet")
all_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset_all.parquet"
)
model_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset.parquet"
)

if (!file.exists(event_file)) {
  stop("Moment-tensor catalogue not found. Run Script 03 first.")
}

if (!file.exists(felt_file)) {
  stop("Felt RAPID observations not found. Run Script 02 first.")
}

events <- read_parquet(event_file) |>
  as_tibble()
felt <- read_parquet(felt_file) |>
  as_tibble() |>
  distinct()

if (!"Mw" %in% names(events)) {
  stop("The Script 03 catalogue does not contain moment-tensor Mw values.")
}

events <- events |>
  rename(
    event_longitude = longitude,
    event_latitude = latitude,
    event_depth_km = depth,
    geonet_catalogue_magnitude = magnitude
  ) |>
  mutate(
    Mw = suppressWarnings(as.numeric(Mw)),
    M_model = Mw,
    M_model_source = if_else(!is.na(Mw), "GeoNet_CMT", NA_character_)
  )

iar_all <- felt |>
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
      !is.na(event_depth_km) & between(event_depth_km, 0, 700)
  )

iar_all$Repi_km <- NA_real_
distance_rows <- which(
  iar_all$valid_report_coordinates & iar_all$valid_event_coordinates
)

if (length(distance_rows) > 0) {
  sf::sf_use_s2(TRUE)

  report_points <- sf::st_as_sf(
    iar_all[distance_rows, ],
    coords = c("report_longitude", "report_latitude"),
    crs = 4326,
    remove = FALSE
  )
  event_points <- sf::st_as_sf(
    iar_all[distance_rows, ],
    coords = c("event_longitude", "event_latitude"),
    crs = 4326,
    remove = FALSE
  )

  iar_all$Repi_km[distance_rows] <- sf::st_distance(
    report_points,
    event_points,
    by_element = TRUE
  ) |>
    units::set_units("km") |>
    units::drop_units() |>
    as.numeric()
}

iar_all <- iar_all |>
  mutate(
    Rhypo_km = if_else(
      !is.na(Repi_km) & valid_event_depth,
      sqrt(Repi_km^2 + event_depth_km^2),
      NA_real_
    ),
    R_model_km = Rhypo_km,
    valid_mmi = !is.na(reported_mmi) & between(reported_mmi, 1, 12),
    valid_report_count =
      !is.na(report_count) & report_count >= minimum_reports_per_cell,
    within_analysis_period =
      !is.na(origintime) &
      origintime >= analysis_start_date &
      (is.na(analysis_end_date) | origintime <= analysis_end_date),
    meets_magnitude_control =
      !is.na(Mw) & Mw > minimum_moment_magnitude,
    within_distance_control =
      !is.na(Repi_km) & Repi_km <= maximum_epicentral_distance_km,
    model_row_eligible =
      valid_mmi &
      valid_report_count &
      valid_report_coordinates &
      valid_event_coordinates &
      valid_event_depth &
      within_analysis_period &
      meets_magnitude_control &
      within_distance_control &
      !is.na(Rhypo_km),
    exclusion_reason = case_when(
      is.na(origintime) ~ "event_not_matched",
      !valid_mmi ~ "invalid_or_missing_mmi",
      !valid_report_count ~ "insufficient_or_missing_report_count",
      !valid_report_coordinates ~ "invalid_report_coordinates",
      !valid_event_coordinates ~ "invalid_event_coordinates",
      !valid_event_depth ~ "invalid_or_missing_event_depth",
      is.na(Mw) ~ "missing_moment_tensor_mw",
      !within_analysis_period ~ "outside_analysis_period",
      !meets_magnitude_control ~ "below_minimum_moment_magnitude",
      !within_distance_control ~ "beyond_maximum_distance",
      is.na(Rhypo_km) ~ "missing_hypocentral_distance",
      TRUE ~ NA_character_
    )
  )

iar <- iar_all |>
  filter(model_row_eligible) |>
  arrange(publicid, Repi_km, cell_id)

if (nrow(iar) == 0) {
  stop("The modelling dataset contains no eligible observations.")
}

dir.create(dirname(model_output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(iar_all, all_output_file)
write_parquet(iar, model_output_file)

message(
  "Saved ", nrow(iar), " modelling rows from ",
  n_distinct(iar$publicid), " earthquakes"
)
