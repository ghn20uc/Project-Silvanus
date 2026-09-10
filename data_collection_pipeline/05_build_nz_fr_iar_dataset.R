# Construct earthquake–Felt RAPID cell dataset

# NOTE TO SELF:
# M_model uses a recognised Mw where available and otherwise uses the original GeoNet catalogue magnitude
# Each Felt cell is compared only with its own earthquake
# Rhypo is calculated using epicentre distance and catalogue hypocentre depth
# Rhypo = sqrt(Repi^2 + depth^2)


library(tidyverse)

# Settings

# Earliest event retained in the modelling-ready dataset
analysis_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)

analysis_end_date <- as.POSIXct(
  NA_character_,
  tz = "UTC"
)

minimum_model_magnitude <- -Inf

maximum_epicentral_distance_km <- Inf

minimum_reports_per_cell <- 1L

minimum_valid_mmi <- 1
maximum_valid_mmi <- 12

# File paths

preferred_event_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt_strong.parquet"
)

fallback_event_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt.parquet"
)

felt_file <- file.path(
  "data_raw",
  "felt_rapid_observations.parquet"
)

all_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset_all.parquet"
)

model_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset.parquet"
)

dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Select event catalogue

if (file.exists(preferred_event_file)) {
  event_catalogue_file <- preferred_event_file
} else if (file.exists(fallback_event_file)) {
  event_catalogue_file <- fallback_event_file
} else {
  stop(
    paste0(
      "No earthquake catalogue was found.\n\n",
      "Expected one of:\n",
      preferred_event_file,
      "\n",
      fallback_event_file,
      "\n\nRun Script 03 and Script 04 first."
    ),
    call. = FALSE
  )
}

if (!file.exists(felt_file)) {
  stop(
    paste0(
      "Felt RAPID observations were not found:\n",
      felt_file,
      "\n\nRun Script 02 first."
    ),
    call. = FALSE
  )
}

# Read inputs

events <- arrow::read_parquet(
  event_catalogue_file
) |>
  as_tibble()
fr <- arrow::read_parquet(
  felt_file
) |>
  as_tibble()
if (nrow(events) == 0) {
  stop(
    "The earthquake catalogue contains zero rows.",
    call. = FALSE
  )
}

if (nrow(fr) == 0) {
  stop(
    "The Felt RAPID observation table contains zero rows.",
    call. = FALSE
  )
}

cat(
  "Earthquakes loaded:",
  format(
    nrow(events),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Felt RAPID cells loaded:",
  format(
    nrow(fr),
    big.mark = ","
  ),
  "\n"
)

# Remove duplicate Felt rows

fr <- fr |>
  distinct()
# Prepare earthquake fields

events <- events |>
  rename(
    event_longitude = longitude,
    event_latitude = latitude,
    event_depth_km = depth,
    catalogue_magnitude = magnitude
  )
# Mw_model should remain a true/explicit Mw.
if (!"Mw_model" %in% names(events)) {
  events$Mw_model <- NA_real_
}

if (!"Mw_source" %in% names(events)) {
  events$Mw_source <- NA_character_
}

events <- events |>
  mutate(
    Mw_model = suppressWarnings(
      as.numeric(Mw_model)
    ),
    M_model = coalesce(
      Mw_model,
      catalogue_magnitude
    ),
    M_model_source = case_when(
      !is.na(Mw_model) &
        !is.na(Mw_source) ~
        Mw_source,
      !is.na(Mw_model) ~
        "recognised_Mw",
      !is.na(catalogue_magnitude) ~
        "GeoNet_catalogue_magnitude",
      TRUE ~
        NA_character_
    ),
    M_model_is_Mw =
      !is.na(Mw_model)
  )

# Join earthquakes data to Felt RAPID cells

iar_all <- fr |>
  left_join(
    events,
    by = "publicid"
  )
# Coordinate validation

iar_all <- iar_all |>
  mutate(
    valid_report_coordinates =
      !is.na(report_longitude) &
      !is.na(report_latitude) &
      report_longitude >= -180 &
      report_longitude <= 180 &
      report_latitude >= -90 &
      report_latitude <= 90,
    valid_event_coordinates =
      !is.na(event_longitude) &
      !is.na(event_latitude) &
      event_longitude >= -180 &
      event_longitude <= 180 &
      event_latitude >= -90 &
      event_latitude <= 90
  )

# Calculate paied epicentral distance

iar_all$Repi_km <- NA_real_

valid_distance_rows <- which(
  iar_all$valid_report_coordinates &
    iar_all$valid_event_coordinates
)

if (length(valid_distance_rows) > 0) {
  sf::sf_use_s2(TRUE)
  report_points <- sf::st_as_sf(
    iar_all[
      valid_distance_rows,
      ,
      drop = FALSE
    ],
    coords = c(
      "report_longitude",
      "report_latitude"
    ),
    crs = 4326,
    remove = FALSE
  )
  event_points <- sf::st_as_sf(
    iar_all[
      valid_distance_rows,
      ,
      drop = FALSE
    ],
    coords = c(
      "event_longitude",
      "event_latitude"
    ),
    crs = 4326,
    remove = FALSE
  )
  paired_distance <- sf::st_distance(
    report_points,
    event_points,
    by_element = TRUE
  )
  iar_all$Repi_km[
    valid_distance_rows
  ] <- paired_distance |>
    units::set_units("km") |>
    units::drop_units() |>
    as.numeric()
}

# Calculate aprox hypocentre distance
iar_all <- iar_all |>
  mutate(
    valid_event_depth =
      !is.na(event_depth_km) &
      event_depth_km >= 0 &
      event_depth_km <= 700,
    Rhypo_km = case_when(
      !is.na(Repi_km) &
        valid_event_depth ~
        sqrt(
          Repi_km^2 +
            event_depth_km^2
        ),
      TRUE ~
        NA_real_
    ),
    Rcentroid_km = NA_real_
  )

if ("mt_centroid_depth_km" %in% names(iar_all)) {
  iar_all <- iar_all |>
    mutate(
      Rcentroid_km = case_when(
        !is.na(Repi_km) &
          !is.na(mt_centroid_depth_km) &
          mt_centroid_depth_km >= 0 &
          mt_centroid_depth_km <= 700 ~
          sqrt(
            Repi_km^2 +
              mt_centroid_depth_km^2
          ),
        TRUE ~
          NA_real_
      )
    )
}

# Default distance 

iar_all <- iar_all |>
  mutate(
    R_model_km = Rhypo_km,
    R_model_type = case_when(
      !is.na(Rhypo_km) ~
        "hypocentral",
      !is.na(Repi_km) ~
        "epicentral_only",
      TRUE ~
        NA_character_
    )
  )
# Modelling elligibility flags 

iar_all <- iar_all |>
  mutate(
    valid_mmi =
      !is.na(reported_mmi) &
      reported_mmi >= minimum_valid_mmi &
      reported_mmi <= maximum_valid_mmi,
    valid_report_count =
      !is.na(report_count) &
      report_count >= minimum_reports_per_cell,
    within_analysis_period =
      !is.na(origintime) &
      origintime >= analysis_start_date &
      (
        is.na(analysis_end_date) |
          origintime <= analysis_end_date
      ),
    meets_magnitude_control =
      !is.na(M_model) &
      M_model >= minimum_model_magnitude,
    within_distance_control =
      !is.na(Repi_km) &
      Repi_km <= maximum_epicentral_distance_km,
    matched_to_event =
      !is.na(origintime),
    model_row_eligible =
      matched_to_event &
      valid_mmi &
      valid_report_count &
      valid_report_coordinates &
      valid_event_coordinates &
      valid_event_depth &
      !is.na(M_model) &
      !is.na(Rhypo_km) &
      within_analysis_period &
      meets_magnitude_control &
      within_distance_control,
    exclusion_reason = case_when(
      !matched_to_event ~
        "event_not_matched",
      !valid_mmi ~
        "invalid_or_missing_mmi",
      !valid_report_count ~
        "insufficient_or_missing_report_count",
      !valid_report_coordinates ~
        "invalid_report_coordinates",
      !valid_event_coordinates ~
        "invalid_event_coordinates",
      !valid_event_depth ~
        "invalid_or_missing_event_depth",
      is.na(M_model) ~
        "missing_model_magnitude",
      is.na(Rhypo_km) ~
        "missing_hypocentral_distance",
      !within_analysis_period ~
        "outside_analysis_period",
      !meets_magnitude_control ~
        "below_minimum_model_magnitude",
      !within_distance_control ~
        "beyond_maximum_distance",
      TRUE ~
        NA_character_
    )
  )
# Create dataset for modelling

iar <- iar_all |>
  filter(
    model_row_eligible
  ) |>
  arrange(
    publicid,
    Repi_km,
    cell_id
  )
if (nrow(iar) == 0) {
  stop(
    paste0(
      "The modelling-ready dataset contains zero rows.\n",
    ),
    call. = FALSE
  )
}

# Save output

# Complete joined dataset
arrow::write_parquet(
  iar_all,
  all_output_file
)

# Modelling-ready dataset
arrow::write_parquet(
  iar,
  model_output_file
)

# Completion summary

cat(
  "\nNZ Felt RAPID attenuation dataset complete\n"
)

