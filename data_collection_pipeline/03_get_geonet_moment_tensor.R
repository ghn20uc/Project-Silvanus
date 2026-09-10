# Download and process the GeoNet moment tensor catalogue, select one preferred solution per earthquake, and join the
# results to the catalogue produced by Script 02.

# NOTES TO SELF:
# moment-tensor data is taken from the GeoNet GitHub https://github.com/GeoNet/data/tree/main/moment-tensor
# Read columns as character first, to prevent old and new ID types conflicting
# Broad fault-style classification based on the direction that one side of a fault moved relative to the other, 
# measured within the fault plane.
# around 0° or ±180°: horizontal movement → strike-slip
# around +90°: upward movement → reverse
# around −90°: downward movement → normal
# intermediate values: mixed → oblique
# the two nodal planes are both saved because the moment tensor alone does not establish which plane is the rupture plane
# a fault style is only assigned if both nodal planes produce the same broad fault-style classification
# Preferred solution ranking is:
# 1. no quality flag
# 2. most complete core CMT information
# 3. highest variance reduction
# 4. greatest station count
# 5. latest solution time
# 6. latest source row

library(tidyverse)

# Set File Paths

felt_catalogue_file <- file.path(
  "data_processed",
  "catalogue_with_felt_summary.parquet"
)

moment_tensor_file <- file.path(
  "data_raw",
  "moment_tensor.parquet"
)

all_solutions_file <- file.path(
  "data_raw",
  "moment_tensor_all_solutions.parquet"
)

enriched_catalogue_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt.parquet"
)

# GeoNet moment tensor source

mt_url <- paste0(
  "https://raw.githubusercontent.com/",
  "GeoNet/data/main/",
  "moment-tensor/GeoNet_CMT_solutions.csv"
)

# Create output directories

dir.create(
  "data_raw",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Check Script 02 output

if (!file.exists(felt_catalogue_file)) {
  stop(
    paste0(
      "The catalogue produced by Script 02 was not found:\n",
      normalizePath(
        felt_catalogue_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 02 before running Script 03."
    ),
    call. = FALSE
  )
}

# Read catalogue containing Felt RAPID summaries

catalogue_with_felt <- arrow::read_parquet(
  felt_catalogue_file
) |>
  as_tibble()
cat(
  "Felt RAPID catalogue loaded:",
  format(
    nrow(catalogue_with_felt),
    big.mark = ","
  ),
  "earthquakes\n"
)

# Download GeoNet moment tensor CSV

download_moment_tensor <- function(url) {
  temporary_file <- tempfile(
    fileext = ".csv"
  )
  on.exit(
    unlink(temporary_file),
    add = TRUE
  )
  cat(
    "Downloading GeoNet moment tensor catalogue...\n"
  )
  response <- tryCatch(
    httr::GET(
      url,
      httr::timeout(60)
    ),
    error = function(e) {
      stop(
        paste0(
          "The GeoNet moment tensor download failed.\n\n",
          "Original error:\n",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  if (httr::http_error(response)) {
    stop(
      paste0(
        "The GeoNet moment tensor download returned HTTP ",
        httr::status_code(response),
        ".\n",
        "URL: ",
        url
      ),
      call. = FALSE
    )
  }
  response_content <- httr::content(
    response,
    as = "raw"
  )
  if (length(response_content) == 0) {
    stop(
      "The downloaded moment tensor file was empty.",
      call. = FALSE
    )
  }
  writeBin(
    response_content,
    temporary_file
  )
  readr::read_csv(
    temporary_file,
    col_types = readr::cols(
      .default = readr::col_character()
    ),
    na = c(
      "",
      "NA",
      "N/A",
      "n/a",
      "NaN",
      "NULL",
      "null",
      "-"
    ),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
}

mt_raw <- download_moment_tensor(
  mt_url
)

if (nrow(mt_raw) == 0) {
  stop(
    "The downloaded moment tensor catalogue contains zero rows.",
    call. = FALSE
  )
}

cat(
  "Moment tensor rows downloaded:",
  format(
    nrow(mt_raw),
    big.mark = ","
  ),
  "\n"
)

# Standardise column names

names(mt_raw) <- names(mt_raw) |>
  stringr::str_trim() |>
  stringr::str_to_lower()
# PublicID or EVENT_ID naming
if (
  "event_id" %in% names(mt_raw) &&
  !"publicid" %in% names(mt_raw)
) {
  mt_raw <- mt_raw |>
    rename(
      publicid = event_id
    )
}

# Check required columns

required_mt_columns <- c(
  "publicid",
  "date",
  "latitude",
  "longitude",
  "strike1",
  "dip1",
  "rake1",
  "strike2",
  "dip2",
  "rake2",
  "ml",
  "mw",
  "mo",
  "cd",
  "ns",
  "dc",
  "vr",
  "method"
)

missing_mt_columns <- setdiff(
  required_mt_columns,
  names(mt_raw)
)

if (length(missing_mt_columns) > 0) {
  stop(
    paste0(
      "The GeoNet moment tensor file is missing required columns:\n",
      paste(
        missing_mt_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        names(mt_raw),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# Parsing

parse_numeric_safe <- function(x) {
  suppressWarnings(
    readr::parse_double(
      as.character(x),
      na = c(
        "",
        "NA",
        "N/A",
        "n/a",
        "NaN",
        "NULL",
        "null",
        "-"
      )
    )
  )
}

parse_mt_time <- function(x) {
  cleaned_time <- x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_replace("\\.0+$", "")
  cleaned_time[
    cleaned_time %in% c(
      "",
      "NA",
      "N/A",
      "NULL",
      "null"
    )
  ] <- NA_character_
  suppressWarnings(
    as.POSIXct(
      cleaned_time,
      format = "%Y%m%d%H%M%S",
      tz = "UTC"
    )
  )
}


classify_rake <- function(rake) {
  normalised_rake <- (
    (rake + 180) %% 360
  ) - 180
  case_when(
    is.na(normalised_rake) ~
      NA_character_,
    abs(normalised_rake) <= 30 ~
      "strike-slip",
    abs(abs(normalised_rake) - 180) <= 30 ~
      "strike-slip",
    normalised_rake > 30 &
      normalised_rake < 150 ~
      "reverse",
    normalised_rake < -30 &
      normalised_rake > -150 ~
      "normal",
    TRUE ~
      "oblique"
  )
}

# Parse numeric moment tensor columns

numeric_mt_columns <- intersect(
  c(
    "latitude",
    "longitude",
    "strike1",
    "dip1",
    "rake1",
    "strike2",
    "dip2",
    "rake2",
    "ml",
    "mw",
    "mo",
    "cd",
    "ns",
    "dc",
    "mxx",
    "mxy",
    "mxz",
    "myy",
    "myz",
    "mzz",
    "vr",
    "tva",
    "tpl",
    "taz",
    "nva",
    "npl",
    "naz",
    "pva",
    "ppl",
    "paz"
  ),
  names(mt_raw)
)

mt <- mt_raw |>
  mutate(
    source_row = row_number(),
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    solution_time = parse_mt_time(
      date
    ),
    across(
      all_of(numeric_mt_columns),
      parse_numeric_safe
    ),
    method = method |>
      as.character() |>
      stringr::str_trim() |>
      na_if("")
  ) |>
  rename(
    mt_latitude = latitude,
    mt_longitude = longitude,
    ML = ml,
    Mw = mw,
    scalar_moment_dyne_cm = mo,
    # GeoNet CMT definitions:
    centroid_depth_km = cd,
    station_count = ns,
    double_couple_percent = dc,
    variance_reduction_percent = vr,
    method_id = method
  ) |>
  mutate(
    station_count = suppressWarnings(
      as.integer(
        round(station_count)
      )
    ),
    # GeoNet uses 9999999 as a placeholder where there is no associated earthquake public ID
    publicid_available =
      !is.na(publicid) &
      publicid != "9999999",
    fault_style_np1 = classify_rake(
      rake1
    ),
    fault_style_np2 = classify_rake(
      rake2
    ),
    fault_style_consensus = case_when(
      is.na(fault_style_np1) |
        is.na(fault_style_np2) ~
        NA_character_,
      fault_style_np1 ==
        fault_style_np2 ~
        fault_style_np1,
      TRUE ~
        NA_character_
    ),
    non_double_couple_percent_derived = case_when(
      !is.na(double_couple_percent) &
        double_couple_percent >= 0 &
        double_couple_percent <= 100 ~
        100 - double_couple_percent,
      TRUE ~
        NA_real_
    ),
    source = "GeoNet_CMT",
    source_url = mt_url,
    downloaded_at_utc = as.POSIXct(
      Sys.time(),
      tz = "UTC"
    )
  ) |>
  select(
    -date
  )
# Data quality flags

mt <- mt |>
  mutate(
    invalid_latitude =
      !is.na(mt_latitude) &
      (
        mt_latitude < -90 |
          mt_latitude > 90
      ),
    invalid_longitude =
      !is.na(mt_longitude) &
      (
        mt_longitude < -180 |
          mt_longitude > 180
      ),
    invalid_centroid_depth =
      !is.na(centroid_depth_km) &
      (
        centroid_depth_km < 0 |
          centroid_depth_km > 700
      ),
    invalid_Mw =
      !is.na(Mw) &
      (
        Mw < 0 |
          Mw > 10
      ),
    invalid_ML =
      !is.na(ML) &
      (
        ML < -2 |
          ML > 10
      ),
    invalid_double_couple =
      !is.na(double_couple_percent) &
      (
        double_couple_percent < 0 |
          double_couple_percent > 100
      ),
    invalid_variance_reduction =
      !is.na(variance_reduction_percent) &
      (
        variance_reduction_percent < -100 |
          variance_reduction_percent > 100
      ),
    invalid_strike =
      (
        !is.na(strike1) &
          (
            strike1 < 0 |
              strike1 > 360
          )
      ) |
      (
        !is.na(strike2) &
          (
            strike2 < 0 |
              strike2 > 360
          )
      ),
    invalid_dip =
      (
        !is.na(dip1) &
          (
            dip1 < 0 |
              dip1 > 90
          )
      ) |
      (
        !is.na(dip2) &
          (
            dip2 < 0 |
              dip2 > 90
          )
      ),
    invalid_rake =
      (
        !is.na(rake1) &
          (
            rake1 < -360 |
              rake1 > 360
          )
      ) |
      (
        !is.na(rake2) &
          (
            rake2 < -360 |
              rake2 > 360
          )
      ),
    quality_flag =
      invalid_latitude |
      invalid_longitude |
      invalid_centroid_depth |
      invalid_Mw |
      invalid_ML |
      invalid_double_couple |
      invalid_variance_reduction |
      invalid_strike |
      invalid_dip |
      invalid_rake
  )
# Choose one preferred solution per earthquake

solution_quality_columns <- intersect(
  c(
    "Mw",
    "scalar_moment_dyne_cm",
    "centroid_depth_km",
    "strike1",
    "dip1",
    "rake1",
    "strike2",
    "dip2",
    "rake2",
    "double_couple_percent",
    "variance_reduction_percent"
  ),
  names(mt)
)

mt_primary <- mt |>
  filter(
    publicid_available
  ) |>
  mutate(
    .completeness_score = rowSums(
      across(
        all_of(solution_quality_columns),
        ~ !is.na(.x)
      )
    ),
    .variance_reduction_rank = coalesce(
      variance_reduction_percent,
      -Inf
    ),
    .station_count_rank = coalesce(
      as.numeric(station_count),
      -Inf
    )
  ) |>
  arrange(
    publicid,
    quality_flag,
    desc(.completeness_score),
    desc(.variance_reduction_rank),
    desc(.station_count_rank),
    desc(solution_time),
    desc(source_row)
  ) |>
  distinct(
    publicid,
    .keep_all = TRUE
  ) |>
  select(
    -.completeness_score,
    -.variance_reduction_rank,
    -.station_count_rank
  )
# Validate primary solution table

if (nrow(mt_primary) == 0) {
  stop(
    paste0(
      "No moment tensor records had a usable publicid.\n",
      "Review ",
      all_solutions_file,
      "."
    ),
    call. = FALSE
  )
}

# Save one-row-per-earthquake CMT table

arrow::write_parquet(
  mt_primary,
  moment_tensor_file
)

# Prepare moment tensor columns for catalogue join

mt_for_join <- mt_primary |>
  transmute(
    publicid,
    has_moment_tensor = TRUE,
    mt_solution_time =
      solution_time,
    mt_latitude,
    mt_longitude,
    mt_ML =
      ML,
    mt_Mw =
      Mw,
    mt_scalar_moment_dyne_cm =
      scalar_moment_dyne_cm,
    mt_centroid_depth_km =
      centroid_depth_km,
    mt_station_count =
      station_count,
    mt_double_couple_percent =
      double_couple_percent,
    mt_non_double_couple_percent_derived =
      non_double_couple_percent_derived,
    mt_variance_reduction_percent =
      variance_reduction_percent,
    mt_strike1 =
      strike1,
    mt_dip1 =
      dip1,
    mt_rake1 =
      rake1,
    mt_strike2 =
      strike2,
    mt_dip2 =
      dip2,
    mt_rake2 =
      rake2,
    mt_fault_style_np1 =
      fault_style_np1,
    mt_fault_style_np2 =
      fault_style_np2,
    mt_fault_style_consensus =
      fault_style_consensus,
    mt_method_id =
      method_id,
    mt_quality_flag =
      quality_flag
  )
# Ensure model columns exist before joining

if (!"Mw_model" %in% names(catalogue_with_felt)) {
  catalogue_with_felt$Mw_model <- NA_real_
}

if (!"fault_style" %in% names(catalogue_with_felt)) {
  catalogue_with_felt$fault_style <- NA_character_
}

# Preserve any values already present in the catalogue
catalogue_with_felt <- catalogue_with_felt |>
  rename(
    Mw_model_before_mt =
      Mw_model,
    fault_style_before_mt =
      fault_style
  )
# Join moment tensor data to Script 02 catalogue

catalogue_with_felt_mt <- catalogue_with_felt |>
  left_join(
    mt_for_join,
    by = "publicid"
  ) |>
  mutate(
    has_moment_tensor = replace_na(
      has_moment_tensor,
      FALSE
    ),
    # Prefer CMT Mw when avilable
    Mw_model = coalesce(
      mt_Mw,
      Mw_model_before_mt
    ),
    Mw_source = case_when(
      !is.na(mt_Mw) ~
        "GeoNet_CMT",
      !is.na(Mw_model_before_mt) ~
        "existing_catalogue_value",
      TRUE ~
        NA_character_
    ),
    fault_style = coalesce(
      mt_fault_style_consensus,
      fault_style_before_mt
    ),
    fault_style_source = case_when(
      !is.na(mt_fault_style_consensus) ~
        "GeoNet_CMT_consensus",
      !is.na(fault_style_before_mt) ~
        "existing_catalogue_value",
      TRUE ~
        NA_character_
    )
  ) |>
  select(
    -Mw_model_before_mt,
    -fault_style_before_mt
  )
# Save  catalogue

arrow::write_parquet(
  catalogue_with_felt_mt,
  enriched_catalogue_file
)

cat(
  "\nMoment tensor processing complete\n"
)

