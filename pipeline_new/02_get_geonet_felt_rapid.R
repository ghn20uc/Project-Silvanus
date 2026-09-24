# Download GeoNet Felt RAPID observations and join query results to the catalogue.

library(tidyverse)
library(arrow)

aggregation_method <- "median"
felt_start_date <- as.POSIXct("2016-09-01 00:00:00", tz = "UTC")
request_pause_seconds <- 0.20
maximum_attempts <- 3L
overwrite_cache <- FALSE

catalogue_file <- file.path("data_raw", "geonet_catalogue.parquet")
cache_dir <- file.path("data_raw", "felt_rapid_cache", aggregation_method)
observation_file <- file.path("data_raw", "felt_rapid_observations.parquet")
catalogue_output_file <- file.path(
  "data_processed",
  "catalogue_with_felt_summary.parquet"
)

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(catalogue_output_file), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(catalogue_file)) {
  stop("GeoNet catalogue not found. Run Script 01 first.")
}

catalogue <- read_parquet(catalogue_file) |>
  as_tibble() |>
  mutate(
    felt_query_eligible =
      !is.na(publicid) &
      !is.na(origintime) &
      origintime >= felt_start_date,
    felt_query_exclusion_reason = case_when(
      is.na(publicid) ~ "missing_publicid",
      is.na(origintime) ~ "missing_origintime",
      origintime < felt_start_date ~ "before_felt_period",
      TRUE ~ NA_character_
    )
  )

query_candidates <- catalogue |>
  filter(felt_query_eligible) |>
  arrange(origintime) |>
  pull(publicid)

if (length(query_candidates) == 0) {
  stop("No earthquakes are eligible for a Felt RAPID query.")
}

empty_observations <- tibble(
  publicid = character(),
  cell_id = character(),
  feature_index = integer(),
  report_longitude = double(),
  report_latitude = double(),
  reported_mmi = double(),
  report_count = integer(),
  geohash = character(),
  aggregation = character()
)

scalar_character <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)

  value <- as.character(unlist(x, recursive = TRUE, use.names = FALSE)[1])
  value <- str_trim(value)

  if (is.na(value) || value == "") NA_character_ else value
}

scalar_numeric <- function(x) {
  value <- scalar_character(x)
  if (is.na(value)) return(NA_real_)
  suppressWarnings(parse_double(value))
}

cache_file_for_id <- function(publicid) {
  safe_id <- str_replace_all(publicid, "[^A-Za-z0-9_-]", "_")
  file.path(cache_dir, paste0(safe_id, ".rds"))
}

parse_feature <- function(feature, publicid, feature_index) {
  properties <- feature$properties
  if (is.null(properties)) properties <- list()
  if (!is.null(names(properties))) names(properties) <- str_to_lower(names(properties))

  coordinates <- unlist(
    feature$geometry$coordinates,
    recursive = TRUE,
    use.names = FALSE
  )

  longitude <- if (length(coordinates) >= 2) as.numeric(coordinates[1]) else NA_real_
  latitude <- if (length(coordinates) >= 2) as.numeric(coordinates[2]) else NA_real_
  geohash <- scalar_character(properties[["geohash"]])
  supplied_id <- scalar_character(feature$id)
  cell_id <- coalesce(supplied_id, geohash, paste0(publicid, "_cell_", feature_index))

  tibble(
    publicid = publicid,
    cell_id = cell_id,
    feature_index = as.integer(feature_index),
    report_longitude = longitude,
    report_latitude = latitude,
    reported_mmi = scalar_numeric(properties[["mmi"]]),
    report_count = as.integer(round(scalar_numeric(properties[["count"]]))),
    geohash = geohash,
    aggregation = aggregation_method
  )
}

make_result <- function(
    publicid,
    query_status,
    http_status = NA_integer_,
    attempt_count = 1L,
    message = NA_character_,
    observations = empty_observations
) {
  list(
    status = tibble(
      publicid = publicid,
      query_status = query_status,
      http_status = as.integer(http_status),
      feature_count = nrow(observations),
      attempt_count = as.integer(attempt_count),
      message = message,
      aggregation = aggregation_method,
      queried_at_utc = as.POSIXct(Sys.time(), tz = "UTC")
    ),
    observations = observations
  )
}

get_felt_rapid <- function(publicid) {
  cache_file <- cache_file_for_id(publicid)

  if (file.exists(cache_file) && !overwrite_cache) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached)) return(cached)
  }

  retryable_status <- c(408L, 429L, 500L, 502L, 503L, 504L)
  last_message <- NA_character_
  last_status <- NA_integer_

  for (attempt in seq_len(maximum_attempts)) {
    Sys.sleep(request_pause_seconds)

    response <- tryCatch(
      httr::GET(
        "https://api.geonet.org.nz/intensity",
        query = list(
          type = "reported",
          aggregation = aggregation_method,
          publicID = publicid
        ),
        httr::add_headers(Accept = "application/vnd.geo+json;version=2"),
        httr::timeout(45)
      ),
      error = function(e) {
        last_message <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(response)) {
      if (attempt < maximum_attempts) {
        Sys.sleep(min(2^attempt, 15))
        next
      }
      result <- make_result(publicid, "network_error", attempt_count = attempt, message = last_message)
      saveRDS(result, cache_file)
      return(result)
    }

    last_status <- httr::status_code(response)

    if (last_status == 204L) {
      result <- make_result(publicid, "no_reports", last_status, attempt)
      saveRDS(result, cache_file)
      return(result)
    }

    if (last_status %in% retryable_status && attempt < maximum_attempts) {
      Sys.sleep(min(2^attempt, 15))
      next
    }

    if (last_status != 200L) {
      response_text <- tryCatch(
        httr::content(response, as = "text", encoding = "UTF-8"),
        error = function(e) NA_character_
      )
      result <- make_result(
        publicid,
        "http_error",
        last_status,
        attempt,
        str_trunc(response_text, 300)
      )
      saveRDS(result, cache_file)
      return(result)
    }

    response_text <- httr::content(response, as = "text", encoding = "UTF-8")
    geojson <- tryCatch(
      jsonlite::fromJSON(response_text, simplifyVector = FALSE),
      error = function(e) {
        last_message <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(geojson)) {
      result <- make_result(
        publicid,
        "json_parse_error",
        last_status,
        attempt,
        last_message
      )
      saveRDS(result, cache_file)
      return(result)
    }

    features <- geojson$features
    if (is.null(features) || length(features) == 0) {
      result <- make_result(publicid, "no_reports", last_status, attempt)
      saveRDS(result, cache_file)
      return(result)
    }

    observations <- map2_dfr(
      features,
      seq_along(features),
      ~ parse_feature(.x, publicid, .y)
    )

    result <- make_result(
      publicid,
      "success",
      last_status,
      attempt,
      observations = observations
    )
    saveRDS(result, cache_file)
    return(result)
  }

  result <- make_result(
    publicid,
    "request_failed",
    last_status,
    maximum_attempts,
    last_message
  )
  saveRDS(result, cache_file)
  result
}

message("Querying Felt RAPID data for ", length(query_candidates), " earthquakes")

for (i in seq_along(query_candidates)) {
  result <- get_felt_rapid(query_candidates[i])

  if (i == 1 || i %% 100 == 0 || i == length(query_candidates)) {
    message(
      "Processed ", i, " of ", length(query_candidates),
      " (", result$status$query_status[1], ")"
    )
  }
}

results <- map(query_candidates, ~ readRDS(cache_file_for_id(.x)))
query_status <- map_dfr(results, "status")
felt_observations <- map_dfr(results, "observations") |>
  distinct()

write_parquet(felt_observations, observation_file)

felt_summary <- felt_observations |>
  group_by(publicid) |>
  summarise(
    felt_cell_count = n(),
    felt_report_count = sum(report_count, na.rm = TRUE),
    .groups = "drop"
  )

query_status <- query_status |>
  transmute(
    publicid,
    felt_query_status = query_status,
    felt_query_http_status = http_status,
    felt_query_attempts = attempt_count,
    felt_query_message = message,
    felt_aggregation = aggregation,
    felt_queried_at_utc = queried_at_utc
  )

catalogue_with_felt <- catalogue |>
  left_join(query_status, by = "publicid") |>
  left_join(felt_summary, by = "publicid") |>
  mutate(
    felt_query_status = if_else(
      felt_query_eligible,
      felt_query_status,
      paste0("not_queried_", felt_query_exclusion_reason)
    ),
    has_felt_rapid = case_when(
      felt_query_status == "success" ~ TRUE,
      felt_query_status == "no_reports" ~ FALSE,
      TRUE ~ NA
    ),
    felt_cell_count = if_else(felt_query_status == "no_reports", 0L, felt_cell_count),
    felt_report_count = if_else(felt_query_status == "no_reports", 0, felt_report_count)
  )

write_parquet(catalogue_with_felt, catalogue_output_file)

message("Saved ", nrow(felt_observations), " Felt RAPID cells to ", observation_file)
