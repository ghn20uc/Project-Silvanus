# ============================================================
# 04_get_geonet_strong_motion.R
#
# Download processed GeoNet strong-motion observations for
# earthquakes in the Felt RAPID and moment-tensor catalogue.
#
# Required input:
#   data_processed/catalogue_with_felt_mt.parquet
#
# Outputs:
#   data_raw/strong_motion_observations.parquet
#   data_raw/strong_motion_query_status.parquet
#   data_processed/catalogue_with_felt_mt_strong.parquet
#
# Cache:
#   data_raw/strong_motion_cache/
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(arrow)


# ============================================================
# DOWNLOAD CONTROLS
# ============================================================

# These settings control which earthquakes are queried from the
# GeoNet API. They are not final model filters.


# Earliest earthquake to query
query_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)


# Latest earthquake to query.
#
# Keep as NA for no upper date limit.
#
# Example:
# query_end_date <- as.POSIXct(
#   "2024-12-31 23:59:59",
#   tz = "UTC"
# )

query_end_date <- as.POSIXct(
  NA_character_,
  tz = "UTC"
)


# Minimum earthquake magnitude to query.
#
# This is an API-query optimisation, not a final model filter.
#
# Suggested workflow:
#   initial run: 3.5
#   expanded run: 3.0
#   maximum coverage: 2.5

minimum_query_magnitude <- 4.0


# TRUE:
#   query only earthquakes that have Felt RAPID observations.
#
# FALSE:
#   query every event satisfying the date and magnitude limits.
#
# TRUE is recommended for this attenuation-relation project,
# because strong-motion observations are primarily needed for
# earthquakes that also have Felt RAPID data.

require_felt_rapid <- TRUE


# Whether events with missing magnitude should still be queried.
#
# FALSE is recommended because these events cannot presently be
# used in a magnitude-dependent attenuation model.

include_missing_magnitude <- FALSE


# Delay between new API requests
request_pause_seconds <- 0.20


# Maximum attempts for temporary API or network failures
maximum_attempts <- 4


# FALSE uses cached results from previous runs.
# TRUE downloads every eligible event again.

overwrite_cache <- FALSE


# Progress-reporting interval
progress_interval <- 100


# ============================================================
# VALIDATE DOWNLOAD CONTROLS
# ============================================================

if (
  length(query_start_date) != 1 ||
  is.na(query_start_date)
) {
  
  stop(
    "query_start_date must contain one valid date-time.",
    call. = FALSE
  )
  
}


if (
  length(query_end_date) != 1
) {
  
  stop(
    "query_end_date must contain one date-time or NA.",
    call. = FALSE
  )
  
}


if (
  !is.na(query_end_date) &&
  query_end_date < query_start_date
) {
  
  stop(
    "query_end_date cannot be earlier than query_start_date.",
    call. = FALSE
  )
  
}


if (
  !is.numeric(minimum_query_magnitude) ||
  length(minimum_query_magnitude) != 1 ||
  is.na(minimum_query_magnitude)
) {
  
  stop(
    "minimum_query_magnitude must be one numeric value.",
    call. = FALSE
  )
  
}


if (
  !is.numeric(request_pause_seconds) ||
  request_pause_seconds < 0
) {
  
  stop(
    "request_pause_seconds must be zero or greater.",
    call. = FALSE
  )
  
}


if (
  !is.numeric(maximum_attempts) ||
  maximum_attempts < 1
) {
  
  stop(
    "maximum_attempts must be at least 1.",
    call. = FALSE
  )
  
}


# ============================================================
# PATHS
# ============================================================

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

status_file <- file.path(
  "data_raw",
  "strong_motion_query_status.parquet"
)

enriched_catalogue_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt_strong.parquet"
)

diagnostic_dir <- file.path(
  "data_processed",
  "diagnostics"
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

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# READ SCRIPT 3 OUTPUT
# ============================================================

if (!file.exists(input_catalogue_file)) {
  
  stop(
    paste0(
      "Script 3 output was not found:\n",
      normalizePath(
        input_catalogue_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 3 before Script 4."
    ),
    call. = FALSE
  )
  
}


catalogue_with_felt_mt <- arrow::read_parquet(
  input_catalogue_file
) |>
  as_tibble()


if (nrow(catalogue_with_felt_mt) == 0) {
  
  stop(
    "The Script 3 catalogue contains zero earthquakes.",
    call. = FALSE
  )
  
}


required_catalogue_columns <- c(
  "publicid",
  "origintime",
  "magnitude"
)


missing_catalogue_columns <- setdiff(
  required_catalogue_columns,
  names(catalogue_with_felt_mt)
)


if (length(missing_catalogue_columns) > 0) {
  
  stop(
    paste0(
      "The Script 3 catalogue is missing required columns:\n",
      paste(
        missing_catalogue_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


if (
  require_felt_rapid &&
  !"has_felt_rapid" %in% names(catalogue_with_felt_mt)
) {
  
  stop(
    paste0(
      "require_felt_rapid is TRUE, but the catalogue does not ",
      "contain a has_felt_rapid column."
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Standardise catalogue fields
# ------------------------------------------------------------

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


catalogue_with_felt_mt <- catalogue_with_felt_mt |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    origintime = parse_origin_time(
      origintime
    ),
    
    magnitude = suppressWarnings(
      as.numeric(magnitude)
    )
    
  )


# Prefer model Mw for query ordering where it is available, but
# fall back to the GeoNet catalogue magnitude.

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


duplicate_catalogue_ids <- catalogue_with_felt_mt |>
  
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


if (nrow(duplicate_catalogue_ids) > 0) {
  
  stop(
    paste0(
      "The input catalogue contains ",
      nrow(duplicate_catalogue_ids),
      " duplicated publicid value(s)."
    ),
    call. = FALSE
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


# ============================================================
# SELECT EARTHQUAKES TO QUERY
# ============================================================

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
      "Review query_start_date, query_end_date, ",
      "minimum_query_magnitude and require_felt_rapid."
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


# ============================================================
# EMPTY TYPED TABLES
# ============================================================

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


empty_status <- tibble(
  
  publicid = character(),
  query_status = character(),
  http_status = integer(),
  
  station_record_count = integer(),
  attempt_count = integer(),
  
  message = character(),
  
  queried_at_utc = as.POSIXct(
    character(),
    tz = "UTC"
  )
  
)


# ============================================================
# VALUE-PARSING HELPERS
# ============================================================

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


# ============================================================
# CACHE HELPERS
# ============================================================

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


cache_is_complete <- function(cache_file) {
  
  if (!file.exists(cache_file)) {
    return(FALSE)
  }
  
  
  cached_result <- read_cache_safely(
    cache_file
  )
  
  
  if (is.null(cached_result)) {
    return(FALSE)
  }
  
  
  cached_status <- cached_result$status$query_status[[1]]
  
  
  cached_status %in% c(
    "success",
    "no_records"
  )
  
}


# ------------------------------------------------------------
# Report cache coverage
# ------------------------------------------------------------

candidate_cache_files <- purrr::map_chr(
  query_candidates$publicid,
  cache_file_for_id
)


complete_cache_count <- sum(
  purrr::map_lgl(
    candidate_cache_files,
    cache_is_complete
  )
)


uncached_count <- if (overwrite_cache) {
  
  nrow(query_candidates)
  
} else {
  
  nrow(query_candidates) -
    complete_cache_count
  
}


minimum_delay_minutes <-
  uncached_count *
  request_pause_seconds /
  60


cat(
  "Completed cached queries:",
  format(
    complete_cache_count,
    big.mark = ","
  ),
  "\n"
)

cat(
  "New or incomplete queries:",
  format(
    uncached_count,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Minimum request-delay time:",
  round(
    minimum_delay_minutes,
    1
  ),
  "minutes\n"
)


# ============================================================
# PARSE ONE STRONG-MOTION FEATURE
# ============================================================

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


# ============================================================
# QUERY-RESULT HELPER
# ============================================================

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


# ============================================================
# DOWNLOAD ONE EARTHQUAKE
# ============================================================

get_strong_motion <- function(publicid) {
  
  cache_file <- cache_file_for_id(
    publicid
  )
  
  
  # Use cached successful/no-record results.
  #
  # Cached failures are retried on later runs.
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
        
        httr::user_agent(
          paste0(
            "DATA309-NZ-Felt-RAPID-IAR/",
            getRversion()
          )
        ),
        
        httr::timeout(45)
        
      ),
      
      error = function(e) {
        
        last_error_message <<-
          conditionMessage(e)
        
        NULL
        
      }
      
    )
    
    
    # --------------------------------------------------------
    # Network-level failure
    # --------------------------------------------------------
    
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
    
    
    # --------------------------------------------------------
    # No strong-motion product
    # --------------------------------------------------------
    
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
    
    
    # --------------------------------------------------------
    # Retry temporary server/API failures
    # --------------------------------------------------------
    
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
    
    
    # --------------------------------------------------------
    # Other HTTP error
    # --------------------------------------------------------
    
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
      
      
      if (
        !is.na(response_text) &&
        nchar(response_text) > 300
      ) {
        
        response_text <- paste0(
          substr(
            response_text,
            1,
            300
          ),
          "..."
        )
        
      }
      
      
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
    
    
    # --------------------------------------------------------
    # Read response body
    # --------------------------------------------------------
    
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
    
    
    # --------------------------------------------------------
    # Parse JSON
    # --------------------------------------------------------
    
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


# ============================================================
# DOWNLOAD ALL ELIGIBLE EARTHQUAKES
# ============================================================

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


# ============================================================
# READ ALL CACHED RESULTS
# ============================================================

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


# ============================================================
# VALIDATE STATION OBSERVATIONS
# ============================================================

if (nrow(strong_motion_observations) > 0) {
  
  invalid_coordinates <- strong_motion_observations |>
    
    filter(
      
      !is.na(station_latitude) &
        (
          station_latitude < -90 |
            station_latitude > 90
        ) |
        
        !is.na(station_longitude) &
        (
          station_longitude < -180 |
            station_longitude > 180
        )
      
    )
  
  
  duplicate_station_records <- strong_motion_observations |>
    
    count(
      publicid,
      station,
      network,
      name = "record_count"
    ) |>
    
    filter(
      record_count > 1
    )
  
} else {
  
  invalid_coordinates <- empty_observations
  
  duplicate_station_records <- tibble(
    publicid = character(),
    station = character(),
    network = character(),
    record_count = integer()
  )
  
}


# ============================================================
# SAVE LONG-FORM STRONG-MOTION DATA
# ============================================================

arrow::write_parquet(
  strong_motion_observations,
  observation_file
)


arrow::write_parquet(
  strong_query_status,
  status_file
)


# ============================================================
# EARTHQUAKE-LEVEL SUMMARY
# ============================================================

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


# ------------------------------------------------------------
# Prepare query status for join
# ------------------------------------------------------------

strong_status_for_join <- strong_query_status |>
  
  transmute(
    
    publicid,
    
    strong_query_status = query_status,
    strong_query_http_status = http_status,
    
    strong_api_station_record_count =
      station_record_count,
    
    strong_query_attempts =
      attempt_count,
    
    strong_query_message =
      message,
    
    strong_queried_at_utc =
      queried_at_utc
    
  )


# ============================================================
# JOIN SUMMARY TO SCRIPT 3 CATALOGUE
# ============================================================

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


# ============================================================
# VALIDATE JOIN
# ============================================================

if (
  nrow(catalogue_with_felt_mt_strong) !=
  nrow(catalogue_with_felt_mt)
) {
  
  stop(
    paste0(
      "The strong-motion join changed the catalogue row count.\n",
      "Before: ",
      nrow(catalogue_with_felt_mt),
      "\nAfter: ",
      nrow(catalogue_with_felt_mt_strong)
    ),
    call. = FALSE
  )
  
}


if (
  anyDuplicated(
    catalogue_with_felt_mt_strong$publicid[
      !is.na(
        catalogue_with_felt_mt_strong$publicid
      )
    ]
  ) > 0
) {
  
  stop(
    paste0(
      "Duplicate publicid values were created by the ",
      "strong-motion join."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# SAVE ENRICHED CATALOGUE
# ============================================================

arrow::write_parquet(
  catalogue_with_felt_mt_strong,
  enriched_catalogue_file
)


# ============================================================
# DIAGNOSTIC OUTPUTS
# ============================================================

strong_query_failures <- strong_query_status |>
  
  filter(
    !query_status %in% c(
      "success",
      "no_records"
    )
  )


strong_status_summary <- strong_query_status |>
  
  count(
    query_status,
    name = "earthquake_count"
  ) |>
  
  arrange(
    desc(earthquake_count)
  )


readr::write_csv(
  strong_query_failures,
  file.path(
    diagnostic_dir,
    "strong_motion_query_failures.csv"
  )
)


readr::write_csv(
  strong_status_summary,
  file.path(
    diagnostic_dir,
    "strong_motion_query_status_summary.csv"
  )
)


readr::write_csv(
  invalid_coordinates,
  file.path(
    diagnostic_dir,
    "strong_motion_invalid_coordinates.csv"
  )
)


readr::write_csv(
  duplicate_station_records,
  file.path(
    diagnostic_dir,
    "strong_motion_duplicate_station_records.csv"
  )
)


# ============================================================
# VERIFY OUTPUTS
# ============================================================

output_files <- c(
  observation_file,
  status_file,
  enriched_catalogue_file
)


invalid_output_files <- output_files[
  
  !file.exists(output_files) |
    is.na(file.info(output_files)$size) |
    file.info(output_files)$size == 0
  
]


if (length(invalid_output_files) > 0) {
  
  stop(
    paste0(
      "One or more output files were not written correctly:\n",
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

events_with_strong_motion <- strong_query_status |>
  
  filter(
    query_status == "success"
  ) |>
  
  summarise(
    n = n_distinct(publicid)
  ) |>
  
  pull(n)


events_without_strong_motion <- strong_query_status |>
  
  filter(
    query_status == "no_records"
  ) |>
  
  summarise(
    n = n_distinct(publicid)
  ) |>
  
  pull(n)


cat(
  "\nStrong-motion processing complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Catalogue earthquakes:",
  format(
    nrow(catalogue_with_felt_mt),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes eligible for query:",
  format(
    nrow(query_candidates),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes with strong-motion data:",
  format(
    events_with_strong_motion,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes with no strong-motion records:",
  format(
    events_without_strong_motion,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Failed or incomplete queries:",
  format(
    nrow(strong_query_failures),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Station-level observations:",
  format(
    nrow(strong_motion_observations),
    big.mark = ","
  ),
  "\n"
)

cat(
  "\nStrong-motion observations:\n",
  normalizePath(
    observation_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nStrong-motion query status:\n",
  normalizePath(
    status_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nEnriched catalogue:\n",
  normalizePath(
    enriched_catalogue_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)