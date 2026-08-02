# ============================================================
# 05_build_nz_fr_iar_dataset.R
#
# Construct the earthquake–Felt RAPID cell dataset used for
# developing a New Zealand intensity attenuation relation.
#
# Preferred event input:
#   data_processed/catalogue_with_felt_mt_strong.parquet
#
# Fallback event input:
#   data_processed/catalogue_with_felt_mt.parquet
#
# Felt RAPID input:
#   data_raw/felt_rapid_observations.parquet
#
# Outputs:
#   data_processed/NZ_FR_IAR_dataset_all.parquet
#   data_processed/NZ_FR_IAR_dataset.parquet
#   data_processed/diagnostics/iar_exclusion_summary.csv
#   data_processed/diagnostics/iar_duplicate_cells.csv
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(arrow)
library(sf)
library(units)


# ============================================================
# MODELLING DATA CONTROLS
# ============================================================

# These controls determine which rows enter the modelling-ready
# dataset. The complete joined dataset is always retained in
# NZ_FR_IAR_dataset_all.parquet.


# Earliest event retained in the modelling-ready dataset
analysis_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)


# Latest event retained.
# Leave as NA for no upper date limit.

analysis_end_date <- as.POSIXct(
  NA_character_,
  tz = "UTC"
)


# Minimum model magnitude.
#
# Use -Inf to retain all events with a usable magnitude.
# A stricter threshold can be applied later after inspecting
# magnitude and distance coverage.

minimum_model_magnitude <- -Inf


# Maximum epicentral distance retained.
#
# Use Inf initially so the distance distribution can be
# inspected before choosing a defensible modelling limit.

maximum_epicentral_distance_km <- Inf


# Minimum number of reports represented by a Felt RAPID cell.
#
# One retains maximum coverage. A later sensitivity analysis
# could test thresholds such as 2, 3, or 5.

minimum_reports_per_cell <- 1L


# Expected broad MMI range used for basic validity checking.
# This is deliberately broad to avoid unnecessary data loss.

minimum_valid_mmi <- 1
maximum_valid_mmi <- 12


# ============================================================
# PATHS
# ============================================================

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

diagnostic_dir <- file.path(
  "data_processed",
  "diagnostics"
)


dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# SELECT EVENT CATALOGUE
# ============================================================

if (file.exists(preferred_event_file)) {
  
  event_catalogue_file <- preferred_event_file
  
} else if (file.exists(fallback_event_file)) {
  
  event_catalogue_file <- fallback_event_file
  
  warning(
    paste0(
      "Script 4 output was not found. Using Script 3 output:\n",
      fallback_event_file,
      "\n\nStrong-motion summaries will not be present, but they ",
      "are not required for calculating the Felt RAPID ",
      "attenuation dataset."
    ),
    call. = FALSE
  )
  
} else {
  
  stop(
    paste0(
      "No enriched earthquake catalogue was found.\n\n",
      "Expected one of:\n",
      preferred_event_file,
      "\n",
      fallback_event_file,
      "\n\nRun Script 3, and preferably Script 4, first."
    ),
    call. = FALSE
  )
  
}


if (!file.exists(felt_file)) {
  
  stop(
    paste0(
      "Felt RAPID observations were not found:\n",
      felt_file,
      "\n\nRun Script 2 first."
    ),
    call. = FALSE
  )
  
}


cat(
  "Event catalogue:",
  event_catalogue_file,
  "\n"
)

cat(
  "Felt RAPID observations:",
  felt_file,
  "\n"
)


# ============================================================
# READ INPUTS
# ============================================================

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


# ============================================================
# VALIDATE REQUIRED COLUMNS
# ============================================================

required_event_columns <- c(
  "publicid",
  "origintime",
  "latitude",
  "longitude",
  "depth",
  "magnitude"
)


required_fr_columns <- c(
  "publicid",
  "report_longitude",
  "report_latitude",
  "reported_mmi",
  "report_count"
)


missing_event_columns <- setdiff(
  required_event_columns,
  names(events)
)


missing_fr_columns <- setdiff(
  required_fr_columns,
  names(fr)
)


if (length(missing_event_columns) > 0) {
  
  stop(
    paste0(
      "The event catalogue is missing required columns:\n",
      paste(
        missing_event_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


if (length(missing_fr_columns) > 0) {
  
  stop(
    paste0(
      "The Felt RAPID table is missing required columns:\n",
      paste(
        missing_fr_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# Create a cell identifier if an older Felt RAPID file lacks one
if (!"cell_id" %in% names(fr)) {
  
  fr <- fr |>
    mutate(
      cell_id = paste0(
        publicid,
        "_cell_",
        row_number()
      )
    )
  
}


# ============================================================
# STANDARDISE INPUT TYPES
# ============================================================

parse_origin_time <- function(x) {
  
  if (inherits(x, "POSIXt")) {
    
    return(
      as.POSIXct(
        x,
        tz = "UTC"
      )
    )
    
  }
  
  
  suppressWarnings(
    lubridate::ymd_hms(
      x,
      quiet = TRUE,
      tz = "UTC"
    )
  )
  
}


events <- events |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    origintime = parse_origin_time(
      origintime
    ),
    
    latitude = suppressWarnings(
      as.numeric(latitude)
    ),
    
    longitude = suppressWarnings(
      as.numeric(longitude)
    ),
    
    depth = suppressWarnings(
      as.numeric(depth)
    ),
    
    magnitude = suppressWarnings(
      as.numeric(magnitude)
    )
    
  )


fr <- fr |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    cell_id = cell_id |>
      as.character() |>
      stringr::str_trim(),
    
    report_longitude = suppressWarnings(
      as.numeric(report_longitude)
    ),
    
    report_latitude = suppressWarnings(
      as.numeric(report_latitude)
    ),
    
    reported_mmi = suppressWarnings(
      as.numeric(reported_mmi)
    ),
    
    report_count = suppressWarnings(
      as.integer(report_count)
    )
    
  )


# ============================================================
# CHECK EARTHQUAKE IDENTIFIERS
# ============================================================

duplicate_event_ids <- events |>
  
  filter(
    !is.na(publicid)
  ) |>
  
  count(
    publicid,
    name = "record_count"
  ) |>
  
  filter(
    record_count > 1
  )


if (nrow(duplicate_event_ids) > 0) {
  
  stop(
    paste0(
      "The earthquake catalogue contains ",
      nrow(duplicate_event_ids),
      " duplicate publicid value(s).\n",
      "The event table must contain one row per earthquake."
    ),
    call. = FALSE
  )
  
}


missing_fr_publicid <- sum(
  is.na(fr$publicid)
)


if (missing_fr_publicid > 0) {
  
  warning(
    paste0(
      missing_fr_publicid,
      " Felt RAPID row(s) have no publicid and cannot be ",
      "matched to an earthquake."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# REMOVE EXACT DUPLICATE FELT ROWS
# ============================================================

fr_rows_before_deduplication <- nrow(fr)


fr <- fr |>
  distinct()


exact_duplicate_rows_removed <-
  fr_rows_before_deduplication -
  nrow(fr)


cat(
  "Exact duplicate Felt rows removed:",
  format(
    exact_duplicate_rows_removed,
    big.mark = ","
  ),
  "\n"
)


# ------------------------------------------------------------
# Identify duplicate event-cell identifiers
# ------------------------------------------------------------

duplicate_felt_cells <- fr |>
  
  filter(
    !is.na(publicid),
    !is.na(cell_id)
  ) |>
  
  count(
    publicid,
    cell_id,
    name = "record_count"
  ) |>
  
  filter(
    record_count > 1
  ) |>
  
  arrange(
    desc(record_count),
    publicid,
    cell_id
  )


if (nrow(duplicate_felt_cells) > 0) {
  
  warning(
    paste0(
      nrow(duplicate_felt_cells),
      " duplicated earthquake–cell identifiers remain after ",
      "removing exact duplicate rows."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# PREPARE EARTHQUAKE FIELDS
# ============================================================

events <- events |>
  
  rename(
    
    event_longitude = longitude,
    event_latitude = latitude,
    
    event_depth_km = depth,
    
    catalogue_magnitude = magnitude
    
  )


# Mw_model should remain a true or explicitly identified Mw.
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
    
    # M_model gives the broadest usable magnitude field.
    #
    # It uses a recognised Mw where available and otherwise
    # falls back to the original GeoNet catalogue magnitude.
    #
    # It is deliberately not called Mw when the fallback source
    # may be ML or another catalogue magnitude type.
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


# ============================================================
# JOIN EARTHQUAKE DATA TO FELT RAPID CELLS
# ============================================================

iar_all <- fr |>
  
  left_join(
    events,
    by = "publicid"
  )


if (nrow(iar_all) != nrow(fr)) {
  
  stop(
    paste0(
      "The earthquake join changed the number of Felt RAPID rows.\n",
      "Before join: ",
      nrow(fr),
      "\nAfter join: ",
      nrow(iar_all)
    ),
    call. = FALSE
  )
  
}


unmatched_event_count <- iar_all |>
  
  filter(
    is.na(origintime)
  ) |>
  
  summarise(
    n = n()
  ) |>
  
  pull(n)


if (unmatched_event_count > 0) {
  
  warning(
    paste0(
      unmatched_event_count,
      " Felt RAPID row(s) did not match an earthquake in the ",
      "event catalogue."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# COORDINATE VALIDATION
# ============================================================

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


# ============================================================
# CALCULATE PAIRED EPICENTRAL DISTANCE
#
# Each Felt cell is compared only with its own earthquake.
# This avoids creation of the incorrect full distance matrix
# produced by the original script.
# ============================================================

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


# ============================================================
# CALCULATE APPROXIMATE HYPOCENTRAL DISTANCE
#
# Rhypo is calculated using epicentral distance and catalogue
# hypocentral depth:
#
#   Rhypo = sqrt(Repi^2 + depth^2)
#
# This is a point-source distance, not rupture distance.
# ============================================================
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
    
    # Initially empty; populated below when CMT centroid depth
    # is present in the dataset.
    Rcentroid_km = NA_real_
    
  )


# Calculate centroid distance separately because the CMT field
# may not be present when Script 3 was skipped or changed.

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


# Default distance for the initial point-source attenuation
# model. A later finite-fault module may add Rrup or Rjb.

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


# ============================================================
# MODELLING ELIGIBILITY FLAGS
# ============================================================

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


# ============================================================
# CREATE MODELLING-READY DATASET
# ============================================================

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
      "Review the modelling controls and exclusion summary."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# DIAGNOSTIC SUMMARIES
# ============================================================

iar_exclusion_summary <- iar_all |>
  
  mutate(
    
    row_status = if_else(
      model_row_eligible,
      "included",
      exclusion_reason
    )
    
  ) |>
  
  count(
    row_status,
    name = "observation_count"
  ) |>
  
  mutate(
    
    percentage = 100 *
      observation_count /
      sum(observation_count)
    
  ) |>
  
  arrange(
    desc(observation_count)
  )


event_coverage_summary <- tibble(
  
  measure = c(
    "Felt RAPID rows loaded",
    "Modelling-ready rows",
    "Unique earthquakes loaded",
    "Unique modelling earthquakes",
    "Rows with recognised Mw",
    "Rows using catalogue magnitude fallback",
    "Rows with moment tensor",
    "Rows with strong motion"
  ),
  
  value = c(
    
    nrow(iar_all),
    
    nrow(iar),
    
    n_distinct(
      iar_all$publicid,
      na.rm = TRUE
    ),
    
    n_distinct(
      iar$publicid,
      na.rm = TRUE
    ),
    
    sum(
      iar_all$M_model_is_Mw,
      na.rm = TRUE
    ),
    
    sum(
      !iar_all$M_model_is_Mw &
        !is.na(iar_all$M_model),
      na.rm = TRUE
    ),
    
    if ("has_moment_tensor" %in% names(iar_all)) {
      sum(
        iar_all$has_moment_tensor %in% TRUE,
        na.rm = TRUE
      )
    } else {
      NA_integer_
    },
    
    if ("has_strong_motion" %in% names(iar_all)) {
      sum(
        iar_all$has_strong_motion %in% TRUE,
        na.rm = TRUE
      )
    } else {
      NA_integer_
    }
    
  )
  
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

# Complete joined dataset, including excluded rows and flags
arrow::write_parquet(
  iar_all,
  all_output_file
)


# Modelling-ready dataset
arrow::write_parquet(
  iar,
  model_output_file
)


readr::write_csv(
  iar_exclusion_summary,
  file.path(
    diagnostic_dir,
    "iar_exclusion_summary.csv"
  )
)


readr::write_csv(
  duplicate_felt_cells,
  file.path(
    diagnostic_dir,
    "iar_duplicate_cells.csv"
  )
)


readr::write_csv(
  event_coverage_summary,
  file.path(
    diagnostic_dir,
    "iar_event_coverage_summary.csv"
  )
)


# ============================================================
# VERIFY OUTPUTS
# ============================================================

output_files <- c(
  all_output_file,
  model_output_file,
  file.path(
    diagnostic_dir,
    "iar_exclusion_summary.csv"
  ),
  file.path(
    diagnostic_dir,
    "iar_duplicate_cells.csv"
  ),
  file.path(
    diagnostic_dir,
    "iar_event_coverage_summary.csv"
  )
)


invalid_output_files <- output_files[
  
  !file.exists(output_files) |
    is.na(file.info(output_files)$size) |
    file.info(output_files)$size == 0
  
]


if (length(invalid_output_files) > 0) {
  
  stop(
    paste0(
      "One or more Script 5 outputs were not written correctly:\n",
      paste(
        invalid_output_files,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
  
}


# ============================================================
# COMPLETION SUMMARY
# ============================================================

cat(
  "\nNZ Felt RAPID attenuation dataset complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Joined Felt RAPID rows:",
  format(
    nrow(iar_all),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Modelling-ready observations:",
  format(
    nrow(iar),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique modelling earthquakes:",
  format(
    n_distinct(
      iar$publicid,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Median reports per retained cell:",
  round(
    median(
      iar$report_count,
      na.rm = TRUE
    ),
    1
  ),
  "\n"
)

cat(
  "Epicentral-distance range:",
  round(
    min(
      iar$Repi_km,
      na.rm = TRUE
    ),
    2
  ),
  "to",
  round(
    max(
      iar$Repi_km,
      na.rm = TRUE
    ),
    2
  ),
  "km\n"
)

cat(
  "Hypocentral-distance range:",
  round(
    min(
      iar$Rhypo_km,
      na.rm = TRUE
    ),
    2
  ),
  "to",
  round(
    max(
      iar$Rhypo_km,
      na.rm = TRUE
    ),
    2
  ),
  "km\n"
)

cat(
  "Model-magnitude range:",
  round(
    min(
      iar$M_model,
      na.rm = TRUE
    ),
    2
  ),
  "to",
  round(
    max(
      iar$M_model,
      na.rm = TRUE
    ),
    2
  ),
  "\n"
)

cat(
  "\nComplete joined output:\n",
  normalizePath(
    all_output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nModelling-ready output:\n",
  normalizePath(
    model_output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)