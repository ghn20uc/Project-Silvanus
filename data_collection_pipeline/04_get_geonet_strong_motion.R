# Adds processed strong-motion information

# NOTES TO SELF:
# source is GeoNet API

library(tidyverse)

# Download Settings

# Earliest earthquake to query
query_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)

# Latest earthquake to query.

query_end_date <- as.POSIXct(
  NA_character_,
  tz = "UTC"
)

# Minimum earthquake magnitude to query.

minimum_query_magnitude <- 5.0

require_felt_rapid <- TRUE

# Whether events with missing magnitude should still be queried

include_missing_magnitude <- FALSE

# Delay between new API requests
request_pause_seconds <- 0.20

# Maximum attempts for temporary API or network failures
maximum_attempts <- 3

# Overwrite the cache

overwrite_cache <- FALSE

# Progress interval
progress_interval <- 100

# Set file path

input_catalogue_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt.parquet"
)

cache_dir <- file.path(
  "data_raw",
  "strong_motion_cache"
)

observation_file <- file.path(
  "data_raw",
  "strong_motion_observations.parquet"
)

enriched_catalogue_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt_strong.parquet"
)

dir.create(
  cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(enriched_catalogue_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# read script 03 output

if (!file.exists(input_catalogue_file)) {
  stop(
    paste0(
      "Script 03 output was not found:\n",
      normalizePath(
        input_catalogue_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 03 before Script 04."
    ),
    call. = FALSE
  )
}

catalogue_with_felt_mt <- arrow::read_parquet(
  input_catalogue_file
) |>
  as_tibble()
# Prefer model Mw for query ordering where available

if ("Mw_model" %in% names(catalogue_with_felt_mt)) {
  catalogue_with_felt_mt <- catalogue_with_felt_mt |>
    mutate(
      query_magnitude = coalesce(
        suppressWarnings(
          as.numeric(Mw_model)
        ),
        magnitude
      )
    )
} else {
  catalogue_with_felt_mt <- catalogue_with_felt_mt |>
    mutate(
      query_magnitude = magnitude
    )
}

cat(
  "Catalogue events loaded:",
  format(
    nrow(catalogue_with_felt_mt),
    big.mark = ","
  ),
  "\n"
)

# Select earthquakes to query

catalogue_with_felt_mt <- catalogue_with_felt_mt |>
  mutate(
    after_start_date =
      !is.na(origintime) &
      origintime >= query_start_date,
    before_end_date =
      is.na(query_end_date) |
      (
        !is.na(origintime) &
          origintime <= query_end_date
      ),
    magnitude_eligible = case_when(
      !is.na(query_magnitude) ~
        query_magnitude >= minimum_query_magnitude,
      include_missing_magnitude ~
        TRUE,
      TRUE ~
        FALSE
    ),
    felt_eligible = case_when(
      !require_felt_rapid ~
        TRUE,
      has_felt_rapid %in% TRUE ~
        TRUE,
      TRUE ~
        FALSE
    ),
    strong_query_eligible =
      !is.na(publicid) &
      !is.na(origintime) &
      after_start_date &
      before_end_date &
      magnitude_eligible &
      felt_eligible,
    strong_query_exclusion_reason = case_when(
      is.na(publicid) ~
        "missing_publicid",
      is.na(origintime) ~
        "missing_origintime",
      !after_start_date ~
        "before_start_date",
      !before_end_date ~
        "after_end_date",
      is.na(query_magnitude) &
        !include_missing_magnitude ~
        "missing_magnitude",
      !is.na(query_magnitude) &
        query_magnitude < minimum_query_magnitude ~
        "below_minimum_magnitude",
      require_felt_rapid &
        !has_felt_rapid %in% TRUE ~
        "no_felt_rapid_data",
      TRUE ~
        NA_character_
    )
  )
query_candidates <- catalogue_with_felt_mt |>
  filter(
    strong_query_eligible
  ) |>
  arrange(
    desc(query_magnitude),
    origintime
  ) |>
  select(
    publicid,
    origintime,
    magnitude,
    query_magnitude,
    has_felt_rapid = any_of("has_felt_rapid")
  )
if (nrow(query_candidates) == 0) {
  stop(
    paste0(
      "No earthquakes satisfy the strong-motion query settings.\n",
      "Review query settings"
    ),
    call. = FALSE
  )
}

cat(
  "Strong-motion query candidates:",
  format(
    nrow(query_candidates),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Query start date:",
  format(
    query_start_date,
    tz = "UTC",
    usetz = TRUE
  ),
  "\n"
)

cat(
  "Query end date:",
  ifelse(
    is.na(query_end_date),
    "No upper limit",
    format(
      query_end_date,
      tz = "UTC",
      usetz = TRUE
    )
  ),
  "\n"
)

cat(
  "Minimum query magnitude:",
  minimum_query_magnitude,
  "\n"
)

cat(
  "Require Felt RAPID data:",
  require_felt_rapid,
  "\n"
)

# Empty typed tables

empty_observations <- tibble(
  publicid = character(),
  feature_index = integer(),
  station = character(),
  station_name = character(),
  network = character(),
  station_location = character(),
  station_longitude = double(),
  station_latitude = double(),
  geometry_type = character(),
  epicentral_distance_km = double(),
  measured_mmi = double(),
  pga_h_percent_g = double(),
  pga_v_percent_g = double(),
  pgv_h_cm_s = double(),
  pgv_v_cm_s = double(),
  strong_event_id = character(),
  strong_event_magnitude = double(),
  strong_event_magnitude_type = character(),
  strong_event_depth_km = double(),
  strong_event_latitude = double(),
  strong_event_longitude = double(),
  strong_product_author = character(),
  strong_product_description = character(),
  strong_product_version = character()
)

# Value parsing

scalar_character <- function(x) {
  if (
    is.null(x) ||
    length(x) == 0
  ) {
    return(NA_character_)
  }
  value <- unlist(
    x,
    recursive = TRUE,
    use.names = FALSE
  )
  if (length(value) == 0) {
    return(NA_character_)
  }
  value <- as.character(
    value[[1]]
  )
  if (
    is.na(value) ||
    stringr::str_trim(value) == ""
  ) {
    return(NA_character_)
  }
  stringr::str_trim(value)
}

scalar_numeric <- function(x) {
  value <- scalar_character(x)
  if (is.na(value)) {
    return(NA_real_)
  }
  suppressWarnings(
    readr::parse_double(
      value,
      na = c(
        "",
        "NA",
        "NaN",
        "NULL",
        "null"
      )
    )
  )
}

# Cache

cache_file_for_id <- function(publicid) {
  safe_id <- stringr::str_replace_all(
    publicid,
    "[^A-Za-z0-9_-]",
    "_"
  )
  file.path(
    cache_dir,
    paste0(
      safe_id,
      ".rds"
    )
  )
}

read_cache_safely <- function(cache_file) {
  tryCatch(
    readRDS(
      cache_file
    ),
    error = function(e) {
      NULL
    }
  )
}

# Parse one strong-motion feature

parse_strong_feature <- function(
    feature,
    publicid,
    feature_index,
    metadata
) {
  if (
    is.null(feature) ||
    !is.list(feature)
  ) {
    return(empty_observations)
  }
  properties <- feature$properties
  if (is.null(properties)) {
    properties <- list()
  }
  coordinates <- unlist(
    feature$geometry$coordinates,
    recursive = TRUE,
    use.names = FALSE
  )
  station_longitude <- NA_real_
  station_latitude <- NA_real_
  if (length(coordinates) >= 2) {
    station_longitude <- suppressWarnings(
      as.numeric(coordinates[[1]])
    )
    station_latitude <- suppressWarnings(
      as.numeric(coordinates[[2]])
    )
  }
  tibble(
    publicid = as.character(publicid),
    feature_index = as.integer(feature_index),
    station = scalar_character(
      properties$station
    ),
    station_name = scalar_character(
      properties$name
    ),
    network = scalar_character(
      properties$network
    ),
    station_location = scalar_character(
      properties$location
    ),
    station_longitude = station_longitude,
    station_latitude = station_latitude,
    geometry_type = scalar_character(
      feature$geometry$type
    ),
    epicentral_distance_km = scalar_numeric(
      properties$distance
    ),
    measured_mmi = scalar_numeric(
      properties$mmi
    ),
    pga_h_percent_g = scalar_numeric(
      properties$pga_h
    ),
    pga_v_percent_g = scalar_numeric(
      properties$pga_v
    ),
    pgv_h_cm_s = scalar_numeric(
      properties$pgv_h
    ),
    pgv_v_cm_s = scalar_numeric(
      properties$pgv_v
    ),
    strong_event_id = scalar_character(
      metadata$event_id
    ),
    strong_event_magnitude = scalar_numeric(
      metadata$magnitude
    ),
    strong_event_magnitude_type = scalar_character(
      metadata$mag_type
    ),
    strong_event_depth_km = scalar_numeric(
      metadata$depth
    ),
    strong_event_latitude = scalar_numeric(
      metadata$latitude
    ),
    strong_event_longitude = scalar_numeric(
      metadata$longitude
    ),
    strong_product_author = scalar_character(
      metadata$author
    ),
    strong_product_description = scalar_character(
      metadata$description
    ),
    strong_product_version = scalar_character(
      metadata$version
    )
  )
}

# Query Result

make_query_result <- function(
    publicid,
    query_status,
    http_status = NA_integer_,
    attempt_count = 1L,
    message = NA_character_,
    observations = empty_observations
) {
  list(
    status = tibble(
      publicid = as.character(publicid),
      query_status = as.character(query_status),
      http_status = as.integer(http_status),
      station_record_count = as.integer(
        nrow(observations)
      ),
      attempt_count = as.integer(
        attempt_count
      ),
      message = as.character(
        message
      ),
      queried_at_utc = as.POSIXct(
        Sys.time(),
        tz = "UTC"
      )
    ),
    observations = observations
  )
}

# Download one earthquake

get_strong_motion <- function(publicid) {
  cache_file <- cache_file_for_id(
    publicid
  )
  if (
    file.exists(cache_file) &&
    !overwrite_cache
  ) {
    cached_result <- read_cache_safely(
      cache_file
    )
    if (!is.null(cached_result)) {
      cached_status <-
        cached_result$status$query_status[[1]]
      if (
        cached_status %in% c(
          "success",
          "no_records"
        )
      ) {
        return(cached_result)
      }
    }
  }
  api_url <- paste0(
    "https://api.geonet.org.nz/",
    "intensity/strong/processed/",
    publicid
  )
  retryable_status_codes <- c(
    408L,
    429L,
    500L,
    502L,
    503L,
    504L
  )
  last_error_message <- NA_character_
  last_http_status <- NA_integer_
  for (attempt in seq_len(maximum_attempts)) {
    Sys.sleep(
      request_pause_seconds
    )
    response <- tryCatch(
      httr::GET(
        api_url,
        httr::add_headers(
          Accept = "application/json",
          `Accept-Encoding` = "gzip"
        ),
        httr::timeout(45)
      ),
      error = function(e) {
        last_error_message <<-
          conditionMessage(e)
        NULL
      }
    )
    if (is.null(response)) {
      if (attempt < maximum_attempts) {
        Sys.sleep(
          min(
            2 ^ attempt,
            15
          )
        )
        next
      }
      result <- make_query_result(
        publicid = publicid,
        query_status = "network_error",
        attempt_count = attempt,
        message = last_error_message
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    http_status <- httr::status_code(
      response
    )
    last_http_status <- http_status
    if (http_status %in% c(204L, 404L)) {
      result <- make_query_result(
        publicid = publicid,
        query_status = "no_records",
        http_status = http_status,
        attempt_count = attempt,
        message = paste0(
          "GeoNet returned HTTP ",
          http_status,
          "."
        )
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    if (
      http_status %in% retryable_status_codes &&
      attempt < maximum_attempts
    ) {
      Sys.sleep(
        min(
          2 ^ attempt,
          15
        )
      )
      next
    }
    if (http_status != 200L) {
      response_text <- tryCatch(
        httr::content(
          response,
          as = "text",
          encoding = "UTF-8"
        ),
        error = function(e) {
          NA_character_
        }
      )
      result <- make_query_result(
        publicid = publicid,
        query_status = "http_error",
        http_status = http_status,
        attempt_count = attempt,
        message = response_text
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    response_text <- tryCatch(
      httr::content(
        response,
        as = "text",
        encoding = "UTF-8"
      ),
      error = function(e) {
        last_error_message <<-
          conditionMessage(e)
        NA_character_
      }
    )
    if (
      is.na(response_text) ||
      stringr::str_trim(response_text) == ""
    ) {
      result <- make_query_result(
        publicid = publicid,
        query_status = "empty_response",
        http_status = http_status,
        attempt_count = attempt,
        message = "GeoNet returned an empty response body."
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    strong_json <- tryCatch(
      jsonlite::fromJSON(
        response_text,
        simplifyVector = FALSE
      ),
      error = function(e) {
        last_error_message <<-
          conditionMessage(e)
        NULL
      }
    )
    if (is.null(strong_json)) {
      result <- make_query_result(
        publicid = publicid,
        query_status = "json_parse_error",
        http_status = http_status,
        attempt_count = attempt,
        message = last_error_message
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    features <- strong_json$features
    metadata <- strong_json$metadata
    if (
      is.null(features) ||
      length(features) == 0
    ) {
      result <- make_query_result(
        publicid = publicid,
        query_status = "no_records",
        http_status = http_status,
        attempt_count = attempt,
        message = "The response contained no station features."
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    observations <- purrr::map2_dfr(
      features,
      seq_along(features),
      function(feature, feature_index) {
        parse_strong_feature(
          feature = feature,
          publicid = publicid,
          feature_index = feature_index,
          metadata = metadata
        )
      }
    )
    if (nrow(observations) == 0) {
      result <- make_query_result(
        publicid = publicid,
        query_status = "feature_parse_error",
        http_status = http_status,
        attempt_count = attempt,
        message = paste0(
          "Station features were returned, but none could be parsed."
        )
      )
      saveRDS(
        result,
        cache_file
      )
      return(result)
    }
    result <- make_query_result(
      publicid = publicid,
      query_status = "success",
      http_status = http_status,
      attempt_count = attempt,
      observations = observations
    )
    saveRDS(
      result,
      cache_file
    )
    return(result)
  }
  result <- make_query_result(
    publicid = publicid,
    query_status = "request_failed",
    http_status = last_http_status,
    attempt_count = maximum_attempts,
    message = last_error_message
  )
  saveRDS(
    result,
    cache_file
  )
  result
}

# Download earthquakes

total_quakes <- nrow(
  query_candidates
)

cat(
  "\nDownloading GeoNet strong-motion data...\n"
)

for (i in seq_len(total_quakes)) {
  current_id <- query_candidates$publicid[[i]]
  result <- get_strong_motion(
    current_id
  )
  if (
    i == 1 ||
    i %% progress_interval == 0 ||
    i == total_quakes
  ) {
    cat(
      "Processed",
      format(i, big.mark = ","),
      "of",
      format(total_quakes, big.mark = ","),
      "|",
      current_id,
      "|",
      result$status$query_status[[1]],
      "\n"
    )
  }
}

cat(
  "Strong-motion download stage complete.\n"
)

# Read all results

read_cached_result <- function(publicid) {
  cache_file <- cache_file_for_id(
    publicid
  )
  cached_result <- read_cache_safely(
    cache_file
  )
  if (!is.null(cached_result)) {
    return(cached_result)
  }
  make_query_result(
    publicid = publicid,
    query_status = "cache_read_error",
    message = "The cached result could not be read."
  )
}

all_results <- purrr::map(
  query_candidates$publicid,
  read_cached_result
)

strong_query_status <- purrr::map_dfr(
  all_results,
  ~ .x$status
)

strong_motion_observations <- purrr::map_dfr(
  all_results,
  ~ .x$observations
)

arrow::write_parquet(
  strong_motion_observations,
  observation_file
)

# Earthquake level summary

max_or_na <- function(x) {
  if (
    length(x) == 0 ||
    all(is.na(x))
  ) {
    return(NA_real_)
  }
  max(
    x,
    na.rm = TRUE
  )
}

median_or_na <- function(x) {
  if (
    length(x) == 0 ||
    all(is.na(x))
  ) {
    return(NA_real_)
  }
  median(
    x,
    na.rm = TRUE
  )
}

if (nrow(strong_motion_observations) > 0) {
  strong_event_summary <- strong_motion_observations |>
    group_by(
      publicid
    ) |>
    summarise(
      strong_station_record_count = n(),
      strong_station_count = n_distinct(
        paste(
          network,
          station,
          sep = "."
        ),
        na.rm = TRUE
      ),
      strong_distance_min_km = if (
        all(is.na(epicentral_distance_km))
      ) {
        NA_real_
      } else {
        min(
          epicentral_distance_km,
          na.rm = TRUE
        )
      },
      strong_distance_max_km = max_or_na(
        epicentral_distance_km
      ),
      strong_mmi_max = max_or_na(
        measured_mmi
      ),
      strong_mmi_median = median_or_na(
        measured_mmi
      ),
      strong_pga_h_max_percent_g = max_or_na(
        pga_h_percent_g
      ),
      strong_pga_v_max_percent_g = max_or_na(
        pga_v_percent_g
      ),
      strong_pgv_h_max_cm_s = max_or_na(
        pgv_h_cm_s
      ),
      strong_pgv_v_max_cm_s = max_or_na(
        pgv_v_cm_s
      ),
      .groups = "drop"
    )
} else {
  strong_event_summary <- tibble(
    publicid = character(),
    strong_station_record_count = integer(),
    strong_station_count = integer(),
    strong_distance_min_km = double(),
    strong_distance_max_km = double(),
    strong_mmi_max = double(),
    strong_mmi_median = double(),
    strong_pga_h_max_percent_g = double(),
    strong_pga_v_max_percent_g = double(),
    strong_pgv_h_max_cm_s = double(),
    strong_pgv_v_max_cm_s = double()
  )
}

# Prepare query status for join

strong_status_for_join <- strong_query_status |>
  transmute(
    publicid,
    strong_query_status = query_status,
    strong_queried_at_utc =
      queried_at_utc
  )
# Join to script 03 output

catalogue_with_felt_mt_strong <-
  catalogue_with_felt_mt |>
  left_join(
    strong_status_for_join,
    by = "publicid"
  ) |>
  left_join(
    strong_event_summary,
    by = "publicid"
  ) |>
  mutate(
    has_strong_motion = case_when(
      strong_query_status == "success" ~
        TRUE,
      strong_query_status == "no_records" ~
        FALSE,
      TRUE ~
        NA
    ),
    strong_query_status = case_when(
      !strong_query_eligible ~
        paste0(
          "not_queried_",
          strong_query_exclusion_reason
        ),
      TRUE ~
        strong_query_status
    ),
    strong_station_record_count = case_when(
      strong_query_status == "no_records" ~
        0L,
      TRUE ~
        strong_station_record_count
    ),
    strong_station_count = case_when(
      strong_query_status == "no_records" ~
        0L,
      TRUE ~
        strong_station_count
    )
  )
# Save catalogue

arrow::write_parquet(
  catalogue_with_felt_mt_strong,
  enriched_catalogue_file
)

cat(
  "\nStrong-motion processing complete\n"
)

