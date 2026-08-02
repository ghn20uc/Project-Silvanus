# ------------------------------------------------------------
# 01_combine_geonet_catalogue.R
#
# Combine GeoNet monthly earthquake catalogue CSV files
# into one validated, unfiltered Parquet catalogue.
# ------------------------------------------------------------

library(tidyverse)
library(arrow)


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_dir <- file.path(
  "data_raw",
  "geonet_catalogue"
)

output_file <- file.path(
  "data_raw",
  "geonet_catalogue.parquet"
)


# ------------------------------------------------------------
# Check input directory
# ------------------------------------------------------------

if (!dir.exists(input_dir)) {
  
  stop(
    paste0(
      "Input directory does not exist:\n",
      normalizePath(
        input_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nCurrent working directory:\n",
      getwd()
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Locate monthly CSV files
# ------------------------------------------------------------

files <- list.files(
  path = input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
) |>
  sort()


if (length(files) == 0) {
  
  stop(
    paste0(
      "No GeoNet monthly CSV files were found in:\n",
      normalizePath(
        input_dir,
        winslash = "/",
        mustWork = TRUE
      )
    ),
    call. = FALSE
  )
  
}


cat(
  "Input directory:",
  normalizePath(
    input_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "Monthly files found:",
  length(files),
  "\n"
)


# ------------------------------------------------------------
# Required GeoNet catalogue columns
# ------------------------------------------------------------

required_columns <- c(
  "publicid",
  "origintime",
  "latitude",
  "longitude",
  "depth",
  "magnitude"
)


# ------------------------------------------------------------
# Read and validate one monthly file
# ------------------------------------------------------------

read_catalogue_file <- function(path) {
  
  file_details <- file.info(path)
  
  if (is.na(file_details$size) || file_details$size == 0) {
    
    stop(
      paste0(
        "Catalogue file is empty or inaccessible:\n",
        path
      ),
      call. = FALSE
    )
    
  }
  
  
  dat <- tryCatch(
    
    readr::read_csv(
      file = path,
      
      # Read everything as character so that differences
      # between monthly files cannot cause bind_rows() errors.
      col_types = readr::cols(
        .default = readr::col_character()
      ),
      
      na = c(
        "",
        "NA",
        "N/A",
        "NaN",
        "NULL",
        "null"
      ),
      
      trim_ws = TRUE,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal"
    ),
    
    error = function(e) {
      
      stop(
        paste0(
          "Failed to read catalogue file:\n",
          path,
          "\n\nOriginal error:\n",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      
    }
    
  )
  
  
  # Standardise capitalisation and remove accidental spaces
  # from column names.
  names(dat) <- names(dat) |>
    stringr::str_trim() |>
    stringr::str_to_lower()
  
  
  duplicated_column_names <- unique(
    names(dat)[duplicated(names(dat))]
  )
  
  
  if (length(duplicated_column_names) > 0) {
    
    stop(
      paste0(
        "Duplicated column names were found in:\n",
        path,
        "\n\nDuplicated columns:\n",
        paste(
          duplicated_column_names,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
    
  }
  
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Required columns are missing from:\n",
        path,
        "\n\nMissing columns:\n",
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
    
  }
  
  
  parsing_problems <- nrow(
    readr::problems(dat)
  )
  
  
  if (parsing_problems > 0) {
    
    warning(
      paste0(
        parsing_problems,
        " malformed CSV row(s) were detected in ",
        basename(path),
        "."
      ),
      call. = FALSE
    )
    
  }
  
  
  dat |>
    mutate(
      source_file = basename(path),
      source_row = row_number(),
      source_mtime = file_details$mtime
    )
  
}


# ------------------------------------------------------------
# Read and combine all monthly files
# ------------------------------------------------------------

raw_catalogue <- purrr::map_dfr(
  files,
  read_catalogue_file
)


if (nrow(raw_catalogue) == 0) {
  
  stop(
    "The monthly files were read, but they contained zero rows.",
    call. = FALSE
  )
  
}


cat(
  "Rows read before cleaning:",
  format(
    nrow(raw_catalogue),
    big.mark = ","
  ),
  "\n"
)


# ------------------------------------------------------------
# Save original values for parsing diagnostics
# ------------------------------------------------------------

raw_origintime <- raw_catalogue$origintime
raw_latitude <- raw_catalogue$latitude
raw_longitude <- raw_catalogue$longitude
raw_depth <- raw_catalogue$depth
raw_magnitude <- raw_catalogue$magnitude


# Helper identifying non-empty source values
has_source_value <- function(x) {
  
  !is.na(x) &
    stringr::str_trim(x) != ""
  
}


# Helper for robust numeric parsing
parse_numeric <- function(x) {
  
  suppressWarnings(
    readr::parse_double(
      x,
      na = c(
        "",
        "NA",
        "N/A",
        "NaN",
        "NULL",
        "null"
      )
    )
  )
  
}


# ------------------------------------------------------------
# Parse important columns
# ------------------------------------------------------------

catalogue_parsed <- raw_catalogue |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    # GeoNet origin times are treated as UTC.
    origintime = suppressWarnings(
      lubridate::ymd_hms(
        origintime,
        tz = "UTC",
        quiet = TRUE
      )
    ),
    
    latitude = parse_numeric(latitude),
    longitude = parse_numeric(longitude),
    depth = parse_numeric(depth),
    magnitude = parse_numeric(magnitude)
    
  )


# ------------------------------------------------------------
# Report failed conversions
# ------------------------------------------------------------

parse_failures <- tibble(
  
  variable = c(
    "origintime",
    "latitude",
    "longitude",
    "depth",
    "magnitude"
  ),
  
  failed_values = c(
    
    sum(
      has_source_value(raw_origintime) &
        is.na(catalogue_parsed$origintime)
    ),
    
    sum(
      has_source_value(raw_latitude) &
        is.na(catalogue_parsed$latitude)
    ),
    
    sum(
      has_source_value(raw_longitude) &
        is.na(catalogue_parsed$longitude)
    ),
    
    sum(
      has_source_value(raw_depth) &
        is.na(catalogue_parsed$depth)
    ),
    
    sum(
      has_source_value(raw_magnitude) &
        is.na(catalogue_parsed$magnitude)
    )
    
  )
  
)


if (any(parse_failures$failed_values > 0)) {
  
  warning(
    paste0(
      "Some non-empty values could not be parsed:\n",
      paste(
        paste0(
          parse_failures$variable,
          ": ",
          parse_failures$failed_values
        ),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Remove rows without an event ID
# ------------------------------------------------------------

missing_publicid_n <- sum(
  is.na(catalogue_parsed$publicid)
)


if (missing_publicid_n > 0) {
  
  warning(
    paste0(
      missing_publicid_n,
      " row(s) had no publicid and were removed."
    ),
    call. = FALSE
  )
  
}


catalogue_parsed <- catalogue_parsed |>
  filter(
    !is.na(publicid)
  )


# ------------------------------------------------------------
# Identify duplicate events
# ------------------------------------------------------------

publicid_counts <- catalogue_parsed |>
  count(
    publicid,
    name = "record_count"
  )


duplicate_event_count <- publicid_counts |>
  filter(
    record_count > 1
  ) |>
  nrow()


duplicate_row_count <- sum(
  publicid_counts$record_count[
    publicid_counts$record_count > 1
  ] - 1
)


if (duplicate_event_count > 0) {
  
  warning(
    paste0(
      duplicate_event_count,
      " publicid value(s) occurred more than once. ",
      duplicate_row_count,
      " duplicate row(s) will be removed."
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Deduplicate
#
# For duplicated public IDs:
#   1. prefer the row with the most complete core data;
#   2. then prefer the most recently modified source file.
# ------------------------------------------------------------

catalogue <- catalogue_parsed |>
  
  mutate(
    
    .core_completeness = rowSums(
      across(
        all_of(
          c(
            "origintime",
            "latitude",
            "longitude",
            "depth",
            "magnitude"
          )
        ),
        ~ !is.na(.x)
      )
    )
    
  ) |>
  
  arrange(
    publicid,
    desc(.core_completeness),
    desc(source_mtime),
    desc(source_row)
  ) |>
  
  distinct(
    publicid,
    .keep_all = TRUE
  ) |>
  
  select(
    -.core_completeness
  )


# ------------------------------------------------------------
# Add derived and placeholder fields
# ------------------------------------------------------------

catalogue <- catalogue |>
  
  mutate(
    
    year = as.integer(
      lubridate::year(origintime)
    ),
    
    # Do not automatically assume catalogue magnitude is Mw.
    Mw_model = NA_real_,
    
    # These should be populated by later enrichment modules.
    tectonic_region = NA_character_,
    fault_style = NA_character_
    
  )


# ------------------------------------------------------------
# Populate Mw_model only where magnitude type explicitly says Mw
# ------------------------------------------------------------

magnitude_type_candidates <- intersect(
  c(
    "magnitudetype",
    "magnitude_type",
    "magtype"
  ),
  names(catalogue)
)


if (length(magnitude_type_candidates) > 0) {
  
  magnitude_type_column <- magnitude_type_candidates[[1]]
  
  magnitude_type_value <- catalogue[[magnitude_type_column]] |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_to_upper()
  
  
  explicit_mw <- magnitude_type_value %in% c(
    "MW",
    "MWB",
    "MWC",
    "MWR",
    "MWW"
  )
  
  
  catalogue$Mw_model[explicit_mw] <-
    catalogue$magnitude[explicit_mw]
  
}


# ------------------------------------------------------------
# Basic coordinate validation
# ------------------------------------------------------------

invalid_latitude_n <- sum(
  !is.na(catalogue$latitude) &
    (
      catalogue$latitude < -90 |
        catalogue$latitude > 90
    )
)


invalid_longitude_n <- sum(
  !is.na(catalogue$longitude) &
    (
      catalogue$longitude < -180 |
        catalogue$longitude > 180
    )
)


if (invalid_latitude_n > 0) {
  
  warning(
    paste0(
      invalid_latitude_n,
      " event(s) have latitude outside -90 to 90."
    ),
    call. = FALSE
  )
  
}


if (invalid_longitude_n > 0) {
  
  warning(
    paste0(
      invalid_longitude_n,
      " event(s) have longitude outside -180 to 180."
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Final catalogue checks
# ------------------------------------------------------------

if (nrow(catalogue) == 0) {
  
  stop(
    "The cleaned catalogue contains zero earthquakes.",
    call. = FALSE
  )
  
}


if (anyDuplicated(catalogue$publicid) > 0) {
  
  stop(
    "Duplicate publicid values remain after deduplication.",
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Write combined catalogue
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


arrow::write_parquet(
  x = catalogue,
  sink = output_file
)


# ------------------------------------------------------------
# Verify output
# ------------------------------------------------------------

if (
  !file.exists(output_file) ||
  is.na(file.info(output_file)$size) ||
  file.info(output_file)$size == 0
) {
  
  stop(
    paste0(
      "The output Parquet file was not created correctly:\n",
      output_file
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------

cat(
  "\nCatalogue successfully combined.\n"
)

cat(
  "Monthly files:",
  length(files),
  "\n"
)

cat(
  "Rows read:",
  format(
    nrow(raw_catalogue),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Rows without publicid removed:",
  format(
    missing_publicid_n,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Duplicate rows removed:",
  format(
    duplicate_row_count,
    big.mark = ","
  ),
  "\n"
)

cat(
  "Final unique events:",
  format(
    nrow(catalogue),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Events with missing magnitude:",
  format(
    sum(is.na(catalogue$magnitude)),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Events with recognised Mw:",
  format(
    sum(!is.na(catalogue$Mw_model)),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Output file:",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "Output size:",
  format(
    file.info(output_file)$size,
    big.mark = ",",
    scientific = FALSE
  ),
  "bytes\n"
)