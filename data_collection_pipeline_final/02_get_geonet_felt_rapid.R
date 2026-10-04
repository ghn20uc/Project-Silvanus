# 02_get_geonet_felt_rapid
# Thomas Rautao 2026

# Download Felt RAPID reports for the catalogue earthquakes

# GeoNet returns reports grouped into geographic cells rather than individual 
# responses. The API is asked for the median MMI in each cell. This
# reduces the influence of unusually high or low individual report while
# retaining the number of reports contributing to each cell value.

library(tidyverse)
library(arrow)

# Median is chosen instead of the default of Max because we want to use a 
# representative cell value rather than the single highest report in each cell.
aggregation_method <- "median"

felt_start_date <- as.POSIXct("2016-09-01 00:00:00", tz = "UTC")

request_pause_seconds <- 0.20
maximum_attempts <- 3L
overwrite_cache <- FALSE

catalogue_file <- file.path("data_raw", "geonet_catalogue.parquet")
cache_dir <- file.path("data_raw", "felt_rapid_cache", aggregation_method)
observation_file <- file.path("data_raw", "felt_rapid_observations.parquet")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(catalogue_file)) {
  stop("GeoNet catalogue not found. Run Script 01 first.")
}

catalogue <- read_parquet(catalogue_file) |>
  as_tibble()

required_columns <- c("publicid", "origintime")
missing_columns <- setdiff(required_columns, names(catalogue))

if (length(missing_columns) > 0) {
  stop("Missing catalogue columns: ", paste(missing_columns, collapse = ", "))
}

query_candidates <- catalogue |>
  filter(
    !is.na(publicid),
    !is.na(origintime),
    origintime >= felt_start_date
  ) |>
  arrange(origintime) |>
  distinct(publicid) |>
  pull(publicid)

if (length(query_candidates) == 0) {
  stop("No earthquakes are eligible for a Felt RAPID query.")
}

# Give every unsuccessful or no-report request the same column structure as 
# a successful request, so they can be combined
empty_observations <- tibble(
  publicid = character(),
  cell_id = character(),
  report_longitude = double(),
  report_latitude = double(),
  reported_mmi = double(),
  report_count = integer()
)

# GeoJSON properties can be missing, nested or represented by different basic
# types. Reduce each property to one clean scalar value.
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

# Rename publicids for safe cache use
cache_file_for_id <- function(publicid) {
  safe_id <- str_replace_all(publicid, "[^A-Za-z0-9_-]", "_")
  file.path(cache_dir, paste0(safe_id, ".rds"))
}

# Convert one GeoJSON feature into the required fields.
# API feature ID is preferred as the cell identifier, then geohash, then 
# a generated identifier
parse_feature <- function(feature, publicid, feature_index) {
  properties <- feature$properties
  if (is.null(properties)) properties <- list()
  if (!is.null(names(properties))) {
    names(properties) <- str_to_lower(names(properties))
  }

  coordinates <- unlist(
    feature$geometry$coordinates,
    recursive = TRUE,
    use.names = FALSE
  )

  longitude <- if (length(coordinates) >= 2) {
    as.numeric(coordinates[1])
  } else {
    NA_real_
  }

  latitude <- if (length(coordinates) >= 2) {
    as.numeric(coordinates[2])
  } else {
    NA_real_
  }

  supplied_id <- scalar_character(feature$id)
  geohash <- scalar_character(properties[["geohash"]])
  cell_id <- coalesce(
    supplied_id,
    geohash,
    paste0(publicid, "_cell_", feature_index)
  )

  tibble(
    publicid = publicid,
    cell_id = cell_id,
    report_longitude = longitude,
    report_latitude = latitude,
    reported_mmi = scalar_numeric(properties[["mmi"]]),
    report_count = as.integer(round(scalar_numeric(properties[["count"]])))
  )
}

make_result <- function(publicid, query_status, observations = empty_observations) {
  list(
    status = tibble(
      publicid = publicid,
      query_status = query_status
    ),
    observations = observations
  )
}

get_felt_rapid <- function(publicid) {
  cache_file <- cache_file_for_id(publicid)

  if (file.exists(cache_file) && !overwrite_cache) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)

    # Only stop for a valid response or a confirmed absence of reports
    cached_status <- tryCatch(
      as.character(cached$status$query_status[1]),
      error = function(e) NA_character_
    )

    if (
      !is.null(cached) &&
        !is.null(cached$observations) &&
        cached_status %in% c("success", "no_reports")
    ) {
      return(cached)
    }
  }

  retryable_status <- c(408L, 429L, 500L, 502L, 503L, 504L)
  last_message <- NA_character_

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

# Failed requests are not cached, so a later run can try them again
    if (is.null(response)) {
      if (attempt < maximum_attempts) {
        Sys.sleep(min(2^attempt, 15))
        next
      }
      return(make_result(publicid, paste0("network_error: ", last_message)))
    }

    http_status <- httr::status_code(response)

    if (http_status == 204L) {
      result <- make_result(publicid, "no_reports")
      saveRDS(result, cache_file)
      return(result)
    }

    if (http_status %in% retryable_status && attempt < maximum_attempts) {
      Sys.sleep(min(2^attempt, 15))
      next
    }

    if (http_status != 200L) {
      return(make_result(publicid, paste0("http_error_", http_status)))
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
      return(make_result(publicid, paste0("json_parse_error: ", last_message)))
    }

    features <- geojson$features
    if (is.null(features) || length(features) == 0) {
      result <- make_result(publicid, "no_reports")
      saveRDS(result, cache_file)
      return(result)
    }

    observations <- map2_dfr(
      features,
      seq_along(features),
      ~ parse_feature(.x, publicid, .y)
    )

    result <- make_result(publicid, "success", observations)
    saveRDS(result, cache_file)
    return(result)
  }

  make_result(publicid, "request_failed")
}

message("Downloading Felt RAPID data for ", length(query_candidates), " earthquakes")

results <- vector("list", length(query_candidates))

for (i in seq_along(query_candidates)) {
  results[[i]] <- get_felt_rapid(query_candidates[i])

  if (i == 1 || i %% 100 == 0 || i == length(query_candidates)) {
    message(
      "Processed ", i, " of ", length(query_candidates),
      " (", results[[i]]$status$query_status[1], ")"
    )
  }
}

query_statuses <- map_chr(
  results,
  ~ as.character(.x$status$query_status[1])
)
incomplete_requests <- !query_statuses %in% c("success", "no_reports")

if (any(incomplete_requests)) {
  incomplete_ids <- query_candidates[incomplete_requests]
  displayed_ids <- head(incomplete_ids, 10)
  remaining_count <- length(incomplete_ids) - length(displayed_ids)
  remaining_text <- if (remaining_count > 0) {
    paste0(" (and ", remaining_count, " more)")
  } else {
    ""
  }

  stop(
    length(incomplete_ids),
    " Felt RAPID requests did not complete. Rerun Script 02; completed ",
    "requests will be read from the cache. Affected IDs: ",
    paste(displayed_ids, collapse = ", "),
    remaining_text
  )
}

felt_observations <- map_dfr(results, "observations") |>
  select(
    publicid,
    cell_id,
    report_longitude,
    report_latitude,
    reported_mmi,
    report_count
  ) |>
  distinct(publicid, cell_id, .keep_all = TRUE) |>
  arrange(publicid, cell_id)

dir.create(dirname(observation_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(felt_observations, observation_file)

message(
  "Saved ", nrow(felt_observations), " Felt RAPID cells from ",
  n_distinct(felt_observations$publicid), " earthquakes"
)
