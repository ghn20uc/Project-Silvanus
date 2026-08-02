# ============================================================
# 02_get_geonet_felt_rapid.R
#
# Download GeoNet Felt RAPID reported-intensity data and attach
# earthquake-level summaries to the GeoNet catalogue.
#
# Inputs:
#   data_raw/geonet_catalogue.parquet
#
# Outputs:
#   data_raw/felt_rapid_observations.parquet
#   data_raw/felt_rapid_mmi_counts.parquet
#   data_raw/felt_rapid_query_status.parquet
#   data_processed/catalogue_with_felt_summary.parquet
#
# Cache:
#   data_raw/felt_rapid_cache/<aggregation>/
# ============================================================


library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(arrow)


# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

# GeoNet supports "median" and "max".
# Median is preferable for a representative intensity within
# each spatial reporting cell.
aggregation_method <- "median"


# Approximate beginning of the Felt RAPID period.
felt_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)


# This controls which earthquakes are queried from the API.
# It is NOT the final model magnitude filter.
#
# Set to -Inf to query every event after felt_start_date.
minimum_query_magnitude <- 4.0


# Delay between uncached requests.
request_pause_seconds <- 0.20


# Maximum attempts for temporary API/network failures.
maximum_attempts <- 4


# Set TRUE only when all cached API results should be replaced.
overwrite_cache <- FALSE


# ------------------------------------------------------------
# Validate settings
# ------------------------------------------------------------

if (!aggregation_method %in% c("median", "max")) {
  
  stop(
    "aggregation_method must be either 'median' or 'max'.",
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


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

catalogue_file <- file.path(
  "data_raw",
  "geonet_catalogue.parquet"
)

cache_dir <- file.path(
  "data_raw",
  "felt_rapid_cache",
  aggregation_method
)

observation_file <- file.path(
  "data_raw",
  "felt_rapid_observations.parquet"
)

mmi_count_file <- file.path(
  "data_raw",
  "felt_rapid_mmi_counts.parquet"
)

status_file <- file.path(
  "data_raw",
  "felt_rapid_query_status.parquet"
)

catalogue_output_file <- file.path(
  "data_processed",
  "catalogue_with_felt_summary.parquet"
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
  dirname(catalogue_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Check catalogue input
# ------------------------------------------------------------

if (!file.exists(catalogue_file)) {
  
  stop(
    paste0(
      "GeoNet catalogue was not found:\n",
      normalizePath(
        catalogue_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 1 first."
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Read catalogue created by Script 1
# ------------------------------------------------------------

catalogue <- arrow::read_parquet(
  catalogue_file
) |>
  as_tibble()


if (nrow(catalogue) == 0) {
  
  stop(
    "The GeoNet catalogue contains zero earthquakes.",
    call. = FALSE
  )
  
}


required_catalogue_columns <- c(
  "publicid",
  "origintime",
  "latitude",
  "longitude",
  "depth",
  "magnitude"
)


missing_catalogue_columns <- setdiff(
  required_catalogue_columns,
  names(catalogue)
)


if (length(missing_catalogue_columns) > 0) {
  
  stop(
    paste0(
      "The GeoNet catalogue is missing required columns:\n",
      paste(
        missing_catalogue_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Standardise catalogue types
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


catalogue <- catalogue |>
  
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


if (anyDuplicated(catalogue$publicid[!is.na(catalogue$publicid)]) > 0) {
  
  stop(
    paste0(
      "Duplicate publicid values remain in the catalogue. ",
      "Script 1 must produce one row per earthquake."
    ),
    call. = FALSE
  )
  
}


cat(
  "Catalogue events loaded:",
  format(
    nrow(catalogue),
    big.mark = ","
  ),
  "\n"
)


# ------------------------------------------------------------
# Identify API-query candidates
#
# Missing-magnitude events are retained because another source
# may later provide their magnitude.
# ------------------------------------------------------------

catalogue <- catalogue |>
  
  mutate(
    
    felt_query_eligible =
      !is.na(publicid) &
      !is.na(origintime) &
      origintime >= felt_start_date &
      (
        is.na(magnitude) |
          magnitude >= minimum_query_magnitude
      ),
    
    felt_query_exclusion_reason = case_when(
      
      is.na(publicid) ~
        "missing_publicid",
      
      is.na(origintime) ~
        "missing_origintime",
      
      origintime < felt_start_date ~
        "before_felt_period",
      
      !is.na(magnitude) &
        magnitude < minimum_query_magnitude ~
        "below_query_magnitude",
      
      TRUE ~
        NA_character_
      
    )
    
  )


query_candidates <- catalogue |>
  
  filter(
    felt_query_eligible
  ) |>
  
  arrange(
    desc(magnitude),
    origintime
  ) |>
  
  select(
    publicid,
    origintime,
    magnitude
  )


if (nrow(query_candidates) == 0) {
  
  stop(
    "No catalogue earthquakes satisfy the Felt RAPID query criteria.",
    call. = FALSE
  )
  
}


cat(
  "Earthquakes eligible for Felt RAPID query:",
  format(
    nrow(query_candidates),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Query period begins:",
  format(
    felt_start_date,
    tz = "UTC",
    usetz = TRUE
  ),
  "\n"
)


cat(
  "Query magnitude threshold:",
  minimum_query_magnitude,
  "\n"
)


cat(
  "Aggregation method:",
  aggregation_method,
  "\n"
)


# ------------------------------------------------------------
# Empty typed outputs
# ------------------------------------------------------------

empty_observations <- tibble(
  
  publicid = character(),
  cell_id = character(),
  feature_index = integer(),
  
  report_longitude = double(),
  report_latitude = double(),
  
  reported_mmi = double(),
  report_count = integer(),
  
  geohash = character(),
  geometry_type = character(),
  
  aggregation = character()
  
)


empty_mmi_counts <- tibble(
  
  publicid = character(),
  cell_id = character(),
  feature_index = integer(),
  
  mmi_level = double(),
  mmi_label = character(),
  report_count_at_mmi = integer(),
  
  aggregation = character()
  
)


empty_status <- tibble(
  
  publicid = character(),
  query_status = character(),
  http_status = integer(),
  
  feature_count = integer(),
  observation_count = integer(),
  mmi_count_rows = integer(),
  
  attempt_count = integer(),
  message = character(),
  
  aggregation = character(),
  queried_at_utc = as.POSIXct(
    character(),
    tz = "UTC"
  )
  
)


# ------------------------------------------------------------
# Scalar conversion helpers
# ------------------------------------------------------------

scalar_character <- function(x) {
  
  if (
    is.null(x) ||
    length(x) == 0
  ) {
    
    return(NA_character_)
    
  }
  
  
  value <- as.character(
    unlist(
      x,
      recursive = TRUE,
      use.names = FALSE
    )[[1]]
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


scalar_integer <- function(x) {
  
  value <- scalar_numeric(x)
  
  
  if (is.na(value)) {
    
    return(NA_integer_)
    
  }
  
  
  as.integer(
    round(value)
  )
  
}


# ------------------------------------------------------------
# Cache filename helper
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Extract MMI-category counts from one feature
# ------------------------------------------------------------

extract_mmi_counts <- function(
    properties,
    publicid,
    cell_id,
    feature_index
) {
  
  if (
    is.null(properties) ||
    length(properties) == 0
  ) {
    
    return(empty_mmi_counts)
    
  }
  
  
  property_names <- names(properties)
  
  
  if (is.null(property_names)) {
    
    return(empty_mmi_counts)
    
  }
  
  
  names(properties) <- property_names |>
    stringr::str_to_lower()
  
  
  count_object <- properties[["count_mmi"]]
  
  
  # Some JSON parsers or API versions may expose these values
  # as separate properties rather than one nested object.
  if (is.null(count_object)) {
    
    matching_names <- names(properties)[
      stringr::str_detect(
        names(properties),
        "^count_mmi[._]"
      )
    ]
    
    
    if (length(matching_names) == 0) {
      
      return(empty_mmi_counts)
      
    }
    
    
    values <- unlist(
      properties[matching_names],
      recursive = TRUE,
      use.names = FALSE
    )
    
    value_names <- matching_names
    
  } else {
    
    values <- unlist(
      count_object,
      recursive = TRUE,
      use.names = TRUE
    )
    
    value_names <- names(values)
    
    
    if (
      is.null(value_names) ||
      all(value_names == "")
    ) {
      
      value_names <- as.character(
        seq_along(values)
      )
      
    }
    
  }
  
  
  if (length(values) == 0) {
    
    return(empty_mmi_counts)
    
  }
  
  
  parsed_counts <- suppressWarnings(
    as.integer(
      as.numeric(values)
    )
  )
  
  
  parsed_levels <- suppressWarnings(
    readr::parse_number(
      value_names
    )
  )
  
  
  result <- tibble(
    
    publicid = publicid,
    cell_id = cell_id,
    feature_index = as.integer(feature_index),
    
    mmi_level = as.numeric(parsed_levels),
    mmi_label = as.character(value_names),
    report_count_at_mmi = parsed_counts,
    
    aggregation = aggregation_method
    
  ) |>
    
    filter(
      !is.na(report_count_at_mmi),
      report_count_at_mmi > 0
    )
  
  
  result
  
}


# ------------------------------------------------------------
# Parse one GeoJSON feature
# ------------------------------------------------------------

parse_intensity_feature <- function(
    feature,
    publicid,
    feature_index
) {
  
  if (
    is.null(feature) ||
    !is.list(feature)
  ) {
    
    return(
      list(
        observation = empty_observations,
        mmi_counts = empty_mmi_counts
      )
    )
    
  }
  
  
  properties <- feature$properties
  
  
  if (is.null(properties)) {
    
    properties <- list()
    
  }
  
  
  if (!is.null(names(properties))) {
    
    names(properties) <- names(properties) |>
      stringr::str_to_lower()
    
  }
  
  
  geometry_type <- scalar_character(
    feature$geometry$type
  )
  
  
  coordinates <- unlist(
    feature$geometry$coordinates,
    recursive = TRUE,
    use.names = FALSE
  )
  
  
  report_longitude <- NA_real_
  report_latitude <- NA_real_
  
  
  if (
    identical(
      geometry_type,
      "Point"
    ) &&
    length(coordinates) >= 2
  ) {
    
    report_longitude <- suppressWarnings(
      as.numeric(coordinates[[1]])
    )
    
    report_latitude <- suppressWarnings(
      as.numeric(coordinates[[2]])
    )
    
  }
  
  
  geohash <- scalar_character(
    properties[["geohash"]]
  )
  
  
  supplied_feature_id <- scalar_character(
    feature$id
  )
  
  
  cell_id <- dplyr::coalesce(
    supplied_feature_id,
    geohash,
    paste0(
      publicid,
      "_cell_",
      feature_index
    )
  )
  
  
  observation <- tibble(
    
    publicid = publicid,
    cell_id = cell_id,
    feature_index = as.integer(feature_index),
    
    report_longitude = report_longitude,
    report_latitude = report_latitude,
    
    reported_mmi = scalar_numeric(
      properties[["mmi"]]
    ),
    
    report_count = scalar_integer(
      properties[["count"]]
    ),
    
    geohash = geohash,
    geometry_type = geometry_type,
    
    aggregation = aggregation_method
    
  )
  
  
  mmi_counts <- extract_mmi_counts(
    properties = properties,
    publicid = publicid,
    cell_id = cell_id,
    feature_index = feature_index
  )
  
  
  list(
    observation = observation,
    mmi_counts = mmi_counts
  )
  
}


# ------------------------------------------------------------
# Create a query-result object
# ------------------------------------------------------------

make_query_result <- function(
    publicid,
    query_status,
    http_status = NA_integer_,
    feature_count = 0L,
    attempt_count = 1L,
    message = NA_character_,
    observations = empty_observations,
    mmi_counts = empty_mmi_counts
) {
  
  list(
    
    status = tibble(
      
      publicid = as.character(publicid),
      query_status = as.character(query_status),
      http_status = as.integer(http_status),
      
      feature_count = as.integer(feature_count),
      observation_count = as.integer(
        nrow(observations)
      ),
      mmi_count_rows = as.integer(
        nrow(mmi_counts)
      ),
      
      attempt_count = as.integer(attempt_count),
      message = as.character(message),
      
      aggregation = aggregation_method,
      queried_at_utc = as.POSIXct(
        Sys.time(),
        tz = "UTC"
      )
      
    ),
    
    observations = observations,
    mmi_counts = mmi_counts
    
  )
  
}


# ------------------------------------------------------------
# Download Felt RAPID data for one earthquake
# ------------------------------------------------------------

get_felt_rapid <- function(publicid) {
  
  cache_file <- cache_file_for_id(
    publicid
  )
  
  
  if (
    file.exists(cache_file) &&
    !overwrite_cache
  ) {
    
    cached_result <- tryCatch(
      
      readRDS(cache_file),
      
      error = function(e) {
        NULL
      }
      
    )
    
    
    if (!is.null(cached_result)) {
      
      return(cached_result)
      
    }
    
  }
  
  
  api_url <- paste0(
    "https://api.geonet.org.nz/intensity"
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
        
        url = api_url,
        
        query = list(
          type = "reported",
          aggregation = aggregation_method,
          publicID = publicid
        ),
        
        httr::add_headers(
          Accept = "application/vnd.geo+json;version=2",
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
        
        last_error_message <<- conditionMessage(e)
        
        NULL
        
      }
      
    )
    
    
    # Network-level failure
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
    
    
    # No-content response
    if (http_status == 204L) {
      
      result <- make_query_result(
        publicid = publicid,
        query_status = "no_reports",
        http_status = http_status,
        attempt_count = attempt,
        message = "GeoNet returned HTTP 204."
      )
      
      
      saveRDS(
        result,
        cache_file
      )
      
      
      return(result)
      
    }
    
    
    # Retry temporary API failures
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
    
    
    # Non-success response
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
    
    
    # Read response body
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
    
    
    # Parse GeoJSON without simplifying nested properties
    geojson <- tryCatch(
      
      jsonlite::fromJSON(
        response_text,
        simplifyVector = FALSE
      ),
      
      error = function(e) {
        
        last_error_message <<- conditionMessage(e)
        
        NULL
        
      }
      
    )
    
    
    if (is.null(geojson)) {
      
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
    
    
    features <- geojson$features
    
    
    if (
      is.null(features) ||
      length(features) == 0
    ) {
      
      result <- make_query_result(
        publicid = publicid,
        query_status = "no_reports",
        http_status = http_status,
        feature_count = 0L,
        attempt_count = attempt,
        message = "GeoJSON contained no intensity features."
      )
      
      
      saveRDS(
        result,
        cache_file
      )
      
      
      return(result)
      
    }
    
    
    parsed_features <- purrr::map2(
      
      features,
      seq_along(features),
      
      function(feature, feature_index) {
        
        parse_intensity_feature(
          feature = feature,
          publicid = publicid,
          feature_index = feature_index
        )
        
      }
      
    )
    
    
    observations <- purrr::map_dfr(
      parsed_features,
      ~ .x$observation
    )
    
    
    mmi_counts <- purrr::map_dfr(
      parsed_features,
      ~ .x$mmi_counts
    )
    
    
    if (nrow(observations) == 0) {
      
      result <- make_query_result(
        publicid = publicid,
        query_status = "feature_parse_error",
        http_status = http_status,
        feature_count = length(features),
        attempt_count = attempt,
        message = paste0(
          "GeoJSON contained features, but none could be parsed."
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
      feature_count = length(features),
      attempt_count = attempt,
      observations = observations,
      mmi_counts = mmi_counts
    )
    
    
    saveRDS(
      result,
      cache_file
    )
    
    
    return(result)
    
  }
  
  
  # This should only be reached if all retry attempts fail.
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


# ------------------------------------------------------------
# Query all eligible earthquakes
#
# Every completed query is written immediately to an RDS cache.
# Therefore the script can be stopped and rerun without starting
# again from the beginning.
# ------------------------------------------------------------

total_quakes <- nrow(
  query_candidates
)


existing_cache_count <- sum(
  file.exists(
    purrr::map_chr(
      query_candidates$publicid,
      cache_file_for_id
    )
  )
)


cat(
  "Previously cached queries:",
  format(
    existing_cache_count,
    big.mark = ","
  ),
  "\n"
)


cat(
  "\nQuerying Felt RAPID data...\n"
)


for (i in seq_len(total_quakes)) {
  
  current_id <- query_candidates$publicid[[i]]
  
  
  result <- get_felt_rapid(
    current_id
  )
  
  
  if (
    i == 1 ||
    i %% 100 == 0 ||
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
  "Felt RAPID querying complete.\n"
)


# ------------------------------------------------------------
# Read cached results
# ------------------------------------------------------------

read_cached_result <- function(publicid) {
  
  cache_file <- cache_file_for_id(
    publicid
  )
  
  
  tryCatch(
    
    readRDS(
      cache_file
    ),
    
    error = function(e) {
      
      make_query_result(
        publicid = publicid,
        query_status = "cache_read_error",
        message = conditionMessage(e)
      )
      
    }
    
  )
  
}


all_results <- purrr::map(
  query_candidates$publicid,
  read_cached_result
)


query_status <- purrr::map_dfr(
  all_results,
  ~ .x$status
)


intensity_observations <- purrr::map_dfr(
  all_results,
  ~ .x$observations
)


mmi_counts <- purrr::map_dfr(
  all_results,
  ~ .x$mmi_counts
)


# ------------------------------------------------------------
# Validate downloaded observations
# ------------------------------------------------------------

if (nrow(intensity_observations) > 0) {
  
  duplicate_cells <- intensity_observations |>
    
    count(
      publicid,
      cell_id,
      name = "record_count"
    ) |>
    
    filter(
      record_count > 1
    )
  
  
  if (nrow(duplicate_cells) > 0) {
    
    warning(
      paste0(
        nrow(duplicate_cells),
        " duplicated event-cell combinations were found."
      ),
      call. = FALSE
    )
    
  }
  
  
  invalid_coordinates <- intensity_observations |>
    
    filter(
      
      !is.na(report_latitude) &
        (
          report_latitude < -90 |
            report_latitude > 90
        ) |
        
        !is.na(report_longitude) &
        (
          report_longitude < -180 |
            report_longitude > 180
        )
      
    )
  
  
  if (nrow(invalid_coordinates) > 0) {
    
    warning(
      paste0(
        nrow(invalid_coordinates),
        " Felt RAPID cells have invalid coordinates."
      ),
      call. = FALSE
    )
    
  }
  
} else {
  
  duplicate_cells <- tibble(
    publicid = character(),
    cell_id = character(),
    record_count = integer()
  )
  
  invalid_coordinates <- empty_observations
  
}


# ------------------------------------------------------------
# Save long-form raw outputs
# ------------------------------------------------------------

arrow::write_parquet(
  intensity_observations,
  observation_file
)


arrow::write_parquet(
  mmi_counts,
  mmi_count_file
)


arrow::write_parquet(
  query_status,
  status_file
)


# ------------------------------------------------------------
# Summary helper functions
# ------------------------------------------------------------

sum_or_na <- function(x) {
  
  if (all(is.na(x))) {
    
    return(NA_real_)
    
  }
  
  
  sum(
    x,
    na.rm = TRUE
  )
  
}


min_or_na <- function(x) {
  
  if (all(is.na(x))) {
    
    return(NA_real_)
    
  }
  
  
  min(
    x,
    na.rm = TRUE
  )
  
}


max_or_na <- function(x) {
  
  if (all(is.na(x))) {
    
    return(NA_real_)
    
  }
  
  
  max(
    x,
    na.rm = TRUE
  )
  
}


median_or_na <- function(x) {
  
  if (all(is.na(x))) {
    
    return(NA_real_)
    
  }
  
  
  median(
    x,
    na.rm = TRUE
  )
  
}


weighted_mean_safe <- function(x, weights) {
  
  valid_x <- is.finite(x)
  
  
  if (!any(valid_x)) {
    
    return(NA_real_)
    
  }
  
  
  valid_weights <- valid_x &
    is.finite(weights) &
    weights > 0
  
  
  if (!any(valid_weights)) {
    
    return(
      mean(
        x[valid_x],
        na.rm = TRUE
      )
    )
    
  }
  
  
  weighted.mean(
    x[valid_weights],
    weights[valid_weights]
  )
  
}


weighted_median_safe <- function(x, weights) {
  
  valid <- is.finite(x) &
    is.finite(weights) &
    weights > 0
  
  
  if (!any(valid)) {
    
    valid_x <- is.finite(x)
    
    
    if (!any(valid_x)) {
      
      return(NA_real_)
      
    }
    
    
    return(
      median(
        x[valid_x]
      )
    )
    
  }
  
  
  x <- x[valid]
  weights <- weights[valid]
  
  
  ordering <- order(x)
  
  x <- x[ordering]
  weights <- weights[ordering]
  
  
  cumulative_weight <- cumsum(weights) /
    sum(weights)
  
  
  x[
    which(
      cumulative_weight >= 0.5
    )[[1]]
  ]
  
}


# ------------------------------------------------------------
# Earthquake-level Felt RAPID summary
# ------------------------------------------------------------

if (nrow(intensity_observations) > 0) {
  
  quake_mmi_summary <- intensity_observations |>
    
    group_by(
      publicid
    ) |>
    
    summarise(
      
      felt_cell_count = n(),
      
      felt_cells_with_mmi = sum(
        !is.na(reported_mmi)
      ),
      
      felt_report_count = sum_or_na(
        report_count
      ),
      
      felt_mmi_min = min_or_na(
        reported_mmi
      ),
      
      felt_mmi_median_cells = median_or_na(
        reported_mmi
      ),
      
      felt_mmi_max = max_or_na(
        reported_mmi
      ),
      
      felt_mmi_weighted_mean = weighted_mean_safe(
        reported_mmi,
        report_count
      ),
      
      felt_mmi_weighted_median = weighted_median_safe(
        reported_mmi,
        report_count
      ),
      
      .groups = "drop"
      
    )
  
} else {
  
  quake_mmi_summary <- tibble(
    
    publicid = character(),
    
    felt_cell_count = integer(),
    felt_cells_with_mmi = integer(),
    felt_report_count = double(),
    
    felt_mmi_min = double(),
    felt_mmi_median_cells = double(),
    felt_mmi_max = double(),
    felt_mmi_weighted_mean = double(),
    felt_mmi_weighted_median = double()
    
  )
  
}


# ------------------------------------------------------------
# Aggregate counts for each MMI category by earthquake
# ------------------------------------------------------------

if (nrow(mmi_counts) > 0) {
  
  quake_mmi_counts <- mmi_counts |>
    
    filter(
      !is.na(mmi_level)
    ) |>
    
    group_by(
      publicid,
      mmi_level
    ) |>
    
    summarise(
      
      report_count_at_mmi = sum(
        report_count_at_mmi,
        na.rm = TRUE
      ),
      
      .groups = "drop"
      
    ) |>
    
    mutate(
      
      mmi_column_level = case_when(
        
        mmi_level == floor(mmi_level) ~
          as.character(
            as.integer(mmi_level)
          ),
        
        TRUE ~
          stringr::str_replace_all(
            as.character(mmi_level),
            "\\.",
            "_"
          )
        
      )
      
    ) |>
    
    select(
      -mmi_level
    ) |>
    
    pivot_wider(
      
      names_from = mmi_column_level,
      values_from = report_count_at_mmi,
      
      names_glue = "felt_mmi_{mmi_column_level}_reports",
      
      values_fill = 0
      
    )
  
} else {
  
  quake_mmi_counts <- tibble(
    publicid = character()
  )
  
}


# ------------------------------------------------------------
# Prepare query status for catalogue join
# ------------------------------------------------------------

query_status_for_join <- query_status |>
  
  transmute(
    
    publicid,
    
    felt_query_status = query_status,
    felt_query_http_status = http_status,
    
    felt_api_feature_count = feature_count,
    felt_observation_count = observation_count,
    
    felt_query_attempts = attempt_count,
    felt_query_message = message,
    
    felt_aggregation = aggregation,
    felt_queried_at_utc = queried_at_utc
    
  )


# ------------------------------------------------------------
# Join Felt RAPID summary to the complete catalogue
#
# has_felt_rapid meanings:
#   TRUE  = queried and reports were found
#   FALSE = queried successfully and no reports were found
#   NA    = not queried, or a request/parsing error occurred
# ------------------------------------------------------------

catalogue_with_felt <- catalogue |>
  
  left_join(
    query_status_for_join,
    by = "publicid"
  ) |>
  
  left_join(
    quake_mmi_summary,
    by = "publicid"
  ) |>
  
  left_join(
    quake_mmi_counts,
    by = "publicid"
  ) |>
  
  mutate(
    
    has_felt_rapid = case_when(
      
      felt_query_status == "success" ~
        TRUE,
      
      felt_query_status == "no_reports" ~
        FALSE,
      
      TRUE ~
        NA
      
    ),
    
    felt_query_status = case_when(
      
      !felt_query_eligible ~
        paste0(
          "not_queried_",
          felt_query_exclusion_reason
        ),
      
      TRUE ~
        felt_query_status
      
    ),
    
    felt_cell_count = case_when(
      
      felt_query_status == "no_reports" ~
        0L,
      
      TRUE ~
        felt_cell_count
      
    ),
    
    felt_cells_with_mmi = case_when(
      
      felt_query_status == "no_reports" ~
        0L,
      
      TRUE ~
        felt_cells_with_mmi
      
    ),
    
    felt_report_count = case_when(
      
      felt_query_status == "no_reports" ~
        0,
      
      TRUE ~
        felt_report_count
      
    )
    
  )


# Set individual MMI report-count columns to zero only where the
# event was successfully queried.
mmi_report_columns <- names(catalogue_with_felt)[
  
  stringr::str_detect(
    names(catalogue_with_felt),
    "^felt_mmi_[0-9_]+_reports$"
  )
  
]


if (length(mmi_report_columns) > 0) {
  
  catalogue_with_felt <- catalogue_with_felt |>
    
    mutate(
      
      across(
        
        all_of(mmi_report_columns),
        
        ~ case_when(
          
          felt_query_status %in% c(
            "success",
            "no_reports"
          ) ~
            replace_na(
              .x,
              0
            ),
          
          TRUE ~
            .x
          
        )
        
      )
      
    )
  
}


# ------------------------------------------------------------
# Save catalogue with Felt RAPID summary
# ------------------------------------------------------------

arrow::write_parquet(
  catalogue_with_felt,
  catalogue_output_file
)


# ------------------------------------------------------------
# Diagnostic outputs
# ------------------------------------------------------------

query_failures <- query_status |>
  
  filter(
    !query_status %in% c(
      "success",
      "no_reports"
    )
  )


readr::write_csv(
  query_failures,
  file.path(
    diagnostic_dir,
    "felt_rapid_query_failures.csv"
  )
)


readr::write_csv(
  duplicate_cells,
  file.path(
    diagnostic_dir,
    "felt_rapid_duplicate_cells.csv"
  )
)


readr::write_csv(
  invalid_coordinates,
  file.path(
    diagnostic_dir,
    "felt_rapid_invalid_coordinates.csv"
  )
)


query_status_summary <- query_status |>
  
  count(
    query_status,
    name = "earthquake_count"
  ) |>
  
  arrange(
    desc(earthquake_count)
  )


readr::write_csv(
  query_status_summary,
  file.path(
    diagnostic_dir,
    "felt_rapid_query_status_summary.csv"
  )
)


# ------------------------------------------------------------
# Verify output files
# ------------------------------------------------------------

output_files <- c(
  observation_file,
  mmi_count_file,
  status_file,
  catalogue_output_file
)


missing_outputs <- output_files[
  !file.exists(output_files)
]


if (length(missing_outputs) > 0) {
  
  stop(
    paste0(
      "The following output files were not created:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------

events_with_reports <- query_status |>
  
  filter(
    query_status == "success"
  ) |>
  
  summarise(
    n = n_distinct(publicid)
  ) |>
  
  pull(n)


events_without_reports <- query_status |>
  
  filter(
    query_status == "no_reports"
  ) |>
  
  summarise(
    n = n_distinct(publicid)
  ) |>
  
  pull(n)


cat(
  "\nFelt RAPID processing complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Catalogue earthquakes:",
  format(
    nrow(catalogue),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes queried:",
  format(
    nrow(query_status),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes with Felt RAPID data:",
  format(
    events_with_reports,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes queried with no reports:",
  format(
    events_without_reports,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Failed or incomplete queries:",
  format(
    nrow(query_failures),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Spatial intensity cells:",
  format(
    nrow(intensity_observations),
    big.mark = ","
  ),
  "\n"
)

cat(
  "MMI count records:",
  format(
    nrow(mmi_counts),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Reported responses represented:",
  format(
    sum(
      intensity_observations$report_count,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "\nObservation output:\n",
  normalizePath(
    observation_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nMMI-count output:\n",
  normalizePath(
    mmi_count_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nQuery-status output:\n",
  normalizePath(
    status_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nCatalogue summary output:\n",
  normalizePath(
    catalogue_output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)