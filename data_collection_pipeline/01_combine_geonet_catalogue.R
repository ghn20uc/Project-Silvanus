# Combine GeoNet monthly CSV files into one Parquet catalogue.

# NOTES TO SELF
# when reading monthly csv files read everything as character to avoid bind_rows() error
# GeoNet origin times are treated as UTC
# can't assume catalogue magnitude is always Mw

library(tidyverse)
library(arrow)


# Set File Paths

input_dir <- file.path(
  "data_raw",
  "geonet_catalogue"
)

output_file <- file.path(
  "data_raw",
  "geonet_catalogue.parquet"
)
# Check input directory

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


# Locate the monthly CSV files

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
  "Monthly files found:",
  length(files),
  "\n"
)
# Read and validate one monthly file

read_catalogue_file <- function(path) {
  
  file_details <- file.info(path)
  dat <- tryCatch(
    
    readr::read_csv(
      file = path,
      
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
  
  
  # Standardise column names.
  names(dat) <- names(dat) |>
    stringr::str_trim() |>
    stringr::str_to_lower()
  dat |>
    mutate(
      source_mtime = file_details$mtime
    )
  
}


# Read and combine all the monthly files

raw_catalogue <- purrr::map_dfr(
  files,
  read_catalogue_file
)


if (nrow(raw_catalogue) == 0) {
  
  stop(
    "Monthly files contained zero rows.",
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
# Numeric parsing
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


# Parse important columns

catalogue_parsed <- raw_catalogue |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
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
# Remove rows without an event ID
catalogue_parsed <- catalogue_parsed |>
  filter(
    !is.na(publicid)
  )
# Deduplicate public IDs
catalogue <- catalogue_parsed |>
  distinct(
    publicid,
    .keep_all = TRUE
  )


# Add placeholder fields

catalogue <- catalogue |>
  
  mutate(
    
    year = as.integer(
      lubridate::year(origintime)
    ),
    
    Mw_model = NA_real_,
    
    # Used by other scripts
    tectonic_region = NA_character_,
    fault_style = NA_character_
    
  )


# Populate Mw_model only where magnitude type explicitly says Mw

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
# Final catalogue checks

if (nrow(catalogue) == 0) {
  
  stop(
    "The cleaned catalogue contains zero earthquakes.",
    call. = FALSE
  )
  
}
# Write combined catalogue

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


arrow::write_parquet(
  x = catalogue,
  sink = output_file
)
cat(
  "\nCatalogue successfully combined.\n"
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
  "Output file:",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
