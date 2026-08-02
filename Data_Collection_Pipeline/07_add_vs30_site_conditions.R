# ============================================================
# 07_add_vs30_site_conditions.R
#
# Extract Foster et al. New Zealand Vs30 raster values at each
# unique Felt RAPID reporting-cell location.
#
# Required input:
#   data_processed/NZ_FR_IAR_with_Rrup.parquet
#
# Required Vs30 source:
#   One GeoTIFF placed in data_raw/vs30/
#
# Outputs:
#   data_processed/vs30_site_lookup.parquet
#   data_processed/NZ_FR_IAR_with_site_conditions.parquet
#   data_processed/NZ_FR_IAR_final_dataset.parquet
#   data_processed/diagnostics/vs30_*.csv
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(arrow)
library(terra)


# ============================================================
# VS30 EXTRACTION CONTROLS
# ============================================================

# Explicitly specify the GeoTIFF path here, or leave as NA to
# automatically select the only .tif/.tiff file in the folder.
#
# Example:
#
# vs30_input_file <- file.path(
#   "data_raw",
#   "vs30",
#   "NZ_Vs30_model.tif"
# )


vs30_input_file <- file.path(
  "data_raw",
  "vs30",
  "combined_vs30.tif"
)

# Directory searched when vs30_input_file is NA
vs30_search_directory <- file.path(
  "data_raw",
  "vs30"
)


# Raster layer to extract.
#
# Use either:
#   1L                  for the first layer
#   "layer_name"        for a named layer

vs30_layer <- 1L


# Extraction method:
#
# "simple"
#   returns the value of the raster cell containing the point.
#
# "bilinear"
#   interpolates from the four surrounding raster cells.
#
# "simple" is recommended initially because it preserves the
# actual published raster-cell values.

vs30_extraction_method <- "simple"


# Where a point falls in a raster cell with an NA value, search
# for the nearest non-missing raster cell within this radius.
#
# Set to 0 to disable nearest-cell filling.
#
# A modest radius is useful for Felt RAPID cell centroids near
# coastlines, but a large value could assign an inappropriate
# site condition from a distant location.

maximum_nearest_search_distance_m <- 2000


# Plausibility limits.
#
# Values outside these limits are retained as Vs30_raw but are
# not used as the cleaned model value.

minimum_plausible_vs30_m_s <- 50
maximum_plausible_vs30_m_s <- 3000


# The main continuous site term will be centred on 760 m/s:
#
#   ln(Vs30 / 760)

vs30_reference_m_s <- 760


# ============================================================
# VALIDATE CONTROLS
# ============================================================

if (
  !vs30_extraction_method %in%
  c("simple", "bilinear")
) {
  
  stop(
    paste0(
      "vs30_extraction_method must be either ",
      "'simple' or 'bilinear'."
    ),
    call. = FALSE
  )
  
}


if (
  !is.numeric(maximum_nearest_search_distance_m) ||
  length(maximum_nearest_search_distance_m) != 1 ||
  is.na(maximum_nearest_search_distance_m) ||
  maximum_nearest_search_distance_m < 0
) {
  
  stop(
    paste0(
      "maximum_nearest_search_distance_m must be one ",
      "non-negative numeric value."
    ),
    call. = FALSE
  )
  
}


if (
  minimum_plausible_vs30_m_s <= 0 ||
  maximum_plausible_vs30_m_s <=
  minimum_plausible_vs30_m_s
) {
  
  stop(
    "The Vs30 plausibility limits are invalid.",
    call. = FALSE
  )
  
}


if (
  !is.numeric(vs30_reference_m_s) ||
  length(vs30_reference_m_s) != 1 ||
  is.na(vs30_reference_m_s) ||
  vs30_reference_m_s <= 0
) {
  
  stop(
    "vs30_reference_m_s must be a positive numeric value.",
    call. = FALSE
  )
  
}


# ============================================================
# PATHS
# ============================================================

input_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_with_Rrup.parquet"
)

site_lookup_file <- file.path(
  "data_processed",
  "vs30_site_lookup.parquet"
)

all_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_with_site_conditions.parquet"
)

final_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_final_dataset.parquet"
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

dir.create(
  vs30_search_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# FIND VS30 RASTER
# ============================================================

if (
  length(vs30_input_file) != 1 ||
  is.na(vs30_input_file) ||
  stringr::str_trim(vs30_input_file) == ""
) {
  
  available_vs30_files <- list.files(
    
    path = vs30_search_directory,
    
    pattern = "\\.(tif|tiff)$",
    
    full.names = TRUE,
    ignore.case = TRUE
    
  ) |>
    sort()
  
  
  if (length(available_vs30_files) == 0) {
    
    stop(
      paste0(
        "No Vs30 GeoTIFF was found in:\n",
        normalizePath(
          vs30_search_directory,
          winslash = "/",
          mustWork = FALSE
        ),
        "\n\nDownload the Foster New Zealand Vs30 GeoTIFF ",
        "and place it in this directory."
      ),
      call. = FALSE
    )
    
  }
  
  
  if (length(available_vs30_files) > 1) {
    
    stop(
      paste0(
        "More than one Vs30 GeoTIFF was found:\n",
        paste(
          available_vs30_files,
          collapse = "\n"
        ),
        "\n\nSet vs30_input_file explicitly at the top ",
        "of Script 7."
      ),
      call. = FALSE
    )
    
  }
  
  
  vs30_input_file <- available_vs30_files[[1]]
  
}


if (!file.exists(vs30_input_file)) {
  
  stop(
    paste0(
      "The specified Vs30 raster does not exist:\n",
      normalizePath(
        vs30_input_file,
        winslash = "/",
        mustWork = FALSE
      )
    ),
    call. = FALSE
  )
  
}


cat(
  "Vs30 raster:",
  normalizePath(
    vs30_input_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)


# ============================================================
# LOAD SCRIPT 6 DATASET
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Script 6 output was not found:\n",
      normalizePath(
        input_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 6 first."
    ),
    call. = FALSE
  )
  
}


iar_rrup <- arrow::read_parquet(
  input_file
) |>
  as_tibble()


if (nrow(iar_rrup) == 0) {
  
  stop(
    "The Script 6 dataset contains zero rows.",
    call. = FALSE
  )
  
}


cat(
  "Felt RAPID modelling observations loaded:",
  format(
    nrow(iar_rrup),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# VALIDATE REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "publicid",
  "report_longitude",
  "report_latitude",
  "reported_mmi"
)


missing_columns <- setdiff(
  required_columns,
  names(iar_rrup)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "The Script 6 dataset is missing required columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# ============================================================
# STANDARDISE COORDINATE TYPES
# ============================================================

iar_rrup <- iar_rrup |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    report_longitude = suppressWarnings(
      as.numeric(report_longitude)
    ),
    
    report_latitude = suppressWarnings(
      as.numeric(report_latitude)
    ),
    
    valid_site_coordinates =
      !is.na(report_longitude) &
      !is.na(report_latitude) &
      report_longitude >= -180 &
      report_longitude <= 180 &
      report_latitude >= -90 &
      report_latitude <= 90
    
  )


invalid_coordinate_count <- sum(
  !iar_rrup$valid_site_coordinates
)


if (invalid_coordinate_count > 0) {
  
  warning(
    paste0(
      invalid_coordinate_count,
      " observation(s) have invalid or missing site ",
      "coordinates and cannot receive Vs30 values."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# LOAD VS30 RASTER
# ============================================================

vs30_raster_all <- tryCatch(
  
  terra::rast(
    vs30_input_file
  ),
  
  error = function(e) {
    
    stop(
      paste0(
        "The Vs30 GeoTIFF could not be opened.\n\n",
        "Original error:\n",
        conditionMessage(e)
      ),
      call. = FALSE
    )
    
  }
  
)


if (terra::nlyr(vs30_raster_all) == 0) {
  
  stop(
    "The Vs30 raster contains no layers.",
    call. = FALSE
  )
  
}


available_layer_names <- names(
  vs30_raster_all
)


# ------------------------------------------------------------
# Select raster layer
# ------------------------------------------------------------

if (is.character(vs30_layer)) {
  
  if (
    length(vs30_layer) != 1 ||
    !vs30_layer %in% available_layer_names
  ) {
    
    stop(
      paste0(
        "The requested Vs30 layer was not found.\n\n",
        "Available layers:\n",
        paste(
          available_layer_names,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
    
  }
  
  
  selected_layer_index <- match(
    vs30_layer,
    available_layer_names
  )
  
} else {
  
  selected_layer_index <- suppressWarnings(
    as.integer(vs30_layer)
  )
  
  
  if (
    length(selected_layer_index) != 1 ||
    is.na(selected_layer_index) ||
    selected_layer_index < 1 ||
    selected_layer_index >
    terra::nlyr(vs30_raster_all)
  ) {
    
    stop(
      paste0(
        "vs30_layer must identify a valid raster layer from 1 to ",
        terra::nlyr(vs30_raster_all),
        "."
      ),
      call. = FALSE
    )
    
  }
  
}


vs30_raster <- vs30_raster_all[[selected_layer_index]]


names(vs30_raster) <- "Vs30_raw"


# ------------------------------------------------------------
# Validate raster CRS
# ------------------------------------------------------------

vs30_crs <- terra::crs(
  vs30_raster
)


if (
  is.na(vs30_crs) ||
  stringr::str_trim(vs30_crs) == ""
) {
  
  stop(
    paste0(
      "The Vs30 raster has no coordinate reference system. ",
      "A valid CRS is required for spatial extraction."
    ),
    call. = FALSE
  )
  
}


cat(
  "Vs30 layer selected:",
  available_layer_names[
    selected_layer_index
  ],
  "\n"
)

cat(
  "Raster dimensions:",
  terra::nrow(vs30_raster),
  "rows ×",
  terra::ncol(vs30_raster),
  "columns\n"
)

cat(
  "Raster resolution:",
  paste(
    round(
      terra::res(vs30_raster),
      4
    ),
    collapse = " × "
  ),
  "\n"
)


# ============================================================
# BUILD UNIQUE SITE LOOKUP
#
# Many earthquakes can contain the same Felt RAPID reporting
# location. Extract Vs30 only once per unique coordinate pair.
# ============================================================

site_lookup <- iar_rrup |>
  
  filter(
    valid_site_coordinates
  ) |>
  
  distinct(
    report_longitude,
    report_latitude
  ) |>
  
  arrange(
    report_longitude,
    report_latitude
  ) |>
  
  mutate(
    site_lookup_id = row_number()
  )


if (nrow(site_lookup) == 0) {
  
  stop(
    "No valid Felt RAPID site coordinates are available.",
    call. = FALSE
  )
  
}


cat(
  "Unique Felt RAPID site locations:",
  format(
    nrow(site_lookup),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# CREATE SPATIAL SITE POINTS
# ============================================================

site_points_wgs84 <- terra::vect(
  
  site_lookup,
  
  geom = c(
    "report_longitude",
    "report_latitude"
  ),
  
  crs = "EPSG:4326",
  keepgeom = TRUE
  
)


site_points_raster_crs <- tryCatch(
  
  terra::project(
    site_points_wgs84,
    vs30_crs
  ),
  
  error = function(e) {
    
    stop(
      paste0(
        "Felt RAPID points could not be transformed to the ",
        "Vs30 raster CRS.\n\n",
        "Original error:\n",
        conditionMessage(e)
      ),
      call. = FALSE
    )
    
  }
  
)


# ============================================================
# EXACT RASTER-CELL EXTRACTION
# ============================================================

exact_extract <- terra::extract(
  
  x = vs30_raster,
  y = site_points_raster_crs,
  
  method = vs30_extraction_method,
  
  cells = TRUE,
  xy = TRUE,
  ID = TRUE,
  
  search_radius = 0
  
)


if (nrow(exact_extract) != nrow(site_lookup)) {
  
  stop(
    paste0(
      "Exact Vs30 extraction returned an unexpected number ",
      "of rows.\n",
      "Sites: ",
      nrow(site_lookup),
      "\nExtraction rows: ",
      nrow(exact_extract)
    ),
    call. = FALSE
  )
  
}


exact_extract <- exact_extract |>
  
  arrange(ID)


site_lookup <- site_lookup |>
  
  mutate(
    
    Vs30_raw = suppressWarnings(
      as.numeric(
        exact_extract$Vs30_raw
      )
    ),
    
    vs30_raster_cell = suppressWarnings(
      as.numeric(
        exact_extract$cell
      )
    ),
    
    vs30_cell_x = suppressWarnings(
      as.numeric(
        exact_extract$x
      )
    ),
    
    vs30_cell_y = suppressWarnings(
      as.numeric(
        exact_extract$y
      )
    ),
    
    vs30_extraction_source = case_when(
      
      !is.na(Vs30_raw) ~
        "containing_raster_cell",
      
      TRUE ~
        "no_raster_value"
      
    )
    
  )


# ============================================================
# NEAREST NON-MISSING CELL FOR EXACT NA VALUES
# ============================================================

missing_exact_indices <- which(
  is.na(site_lookup$Vs30_raw)
)


if (
  length(missing_exact_indices) > 0 &&
  maximum_nearest_search_distance_m > 0
) {
  
  cat(
    "Sites requiring nearest-cell search:",
    format(
      length(missing_exact_indices),
      big.mark = ","
    ),
    "\n"
  )
  
  
  missing_points <- site_points_raster_crs[
    missing_exact_indices
  ]
  
  
  nearest_extract <- terra::extract(
    
    x = vs30_raster,
    y = missing_points,
    
    method = "simple",
    
    cells = TRUE,
    xy = TRUE,
    ID = TRUE,
    
    search_radius =
      maximum_nearest_search_distance_m
    
  ) |>
    as.data.frame() |>
    arrange(ID)
  
  
  if (
    nrow(nearest_extract) !=
    length(missing_exact_indices)
  ) {
    
    stop(
      paste0(
        "Nearest-cell extraction returned an unexpected ",
        "number of rows.\n",
        "Expected: ",
        length(missing_exact_indices),
        "\nReturned: ",
        nrow(nearest_extract)
      ),
      call. = FALSE
    )
    
  }
  
  
  cat(
    "Nearest extraction columns:",
    paste(
      names(nearest_extract),
      collapse = ", "
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Identify returned metadata columns robustly
  # ----------------------------------------------------------
  
  find_first_column <- function(
    column_names,
    patterns
  ) {
    
    for (pattern in patterns) {
      
      matching_columns <- grep(
        pattern,
        column_names,
        value = TRUE,
        ignore.case = TRUE
      )
      
      
      if (length(matching_columns) > 0) {
        
        return(
          matching_columns[[1]]
        )
        
      }
      
    }
    
    
    NA_character_
    
  }
  
  
  cell_column <- find_first_column(
    
    names(nearest_extract),
    
    c(
      "^cell$",
      "^cells$",
      "^cell[._]",
      "^cell_number$",
      "^cellnumber$"
    )
    
  )
  
  
  x_column <- find_first_column(
    
    names(nearest_extract),
    
    c(
      "^x$",
      "^x[._]",
      "^x_coord",
      "^xcoord"
    )
    
  )
  
  
  y_column <- find_first_column(
    
    names(nearest_extract),
    
    c(
      "^y$",
      "^y[._]",
      "^y_coord",
      "^ycoord"
    )
    
  )
  
  
  distance_column <- find_first_column(
    
    names(nearest_extract),
    
    c(
      "^distance$",
      "^distance[._]",
      "^dist$",
      "^dist[._]"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Extract the returned Vs30 values
  # ----------------------------------------------------------
  
  nearest_values <- suppressWarnings(
    as.numeric(
      nearest_extract$Vs30_raw
    )
  )
  
  
  nearest_value_available <-
    !is.na(nearest_values) &
    is.finite(nearest_values)
  
  
  rows_to_update <- missing_exact_indices[
    nearest_value_available
  ]
  
  
  nearest_source_rows <- which(
    nearest_value_available
  )
  
  
  # ----------------------------------------------------------
  # Extract cell numbers when supplied
  # ----------------------------------------------------------
  
  nearest_cells <- rep(
    NA_real_,
    nrow(nearest_extract)
  )
  
  
  if (!is.na(cell_column)) {
    
    nearest_cells <- suppressWarnings(
      as.numeric(
        nearest_extract[[cell_column]]
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Extract cell-centre coordinates when supplied
  # ----------------------------------------------------------
  
  nearest_x <- rep(
    NA_real_,
    nrow(nearest_extract)
  )
  
  
  nearest_y <- rep(
    NA_real_,
    nrow(nearest_extract)
  )
  
  
  if (!is.na(x_column)) {
    
    nearest_x <- suppressWarnings(
      as.numeric(
        nearest_extract[[x_column]]
      )
    )
    
  }
  
  
  if (!is.na(y_column)) {
    
    nearest_y <- suppressWarnings(
      as.numeric(
        nearest_extract[[y_column]]
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Recover cell numbers from returned x/y where necessary
  # ----------------------------------------------------------
  
  cells_missing_but_xy_available <-
    is.na(nearest_cells) &
    !is.na(nearest_x) &
    !is.na(nearest_y)
  
  
  if (any(cells_missing_but_xy_available)) {
    
    nearest_cells[
      cells_missing_but_xy_available
    ] <- terra::cellFromXY(
      
      vs30_raster,
      
      cbind(
        
        nearest_x[
          cells_missing_but_xy_available
        ],
        
        nearest_y[
          cells_missing_but_xy_available
        ]
        
      )
      
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Recover cell-centre x/y from cell numbers where necessary
  # ----------------------------------------------------------
  
  xy_missing_but_cell_available <-
    !is.na(nearest_cells) &
    (
      is.na(nearest_x) |
        is.na(nearest_y)
    )
  
  
  if (any(xy_missing_but_cell_available)) {
    
    recovered_xy <- terra::xyFromCell(
      
      vs30_raster,
      
      nearest_cells[
        xy_missing_but_cell_available
      ]
      
    )
    
    
    nearest_x[
      xy_missing_but_cell_available
    ] <- recovered_xy[, 1]
    
    
    nearest_y[
      xy_missing_but_cell_available
    ] <- recovered_xy[, 2]
    
  }
  
  
  # ----------------------------------------------------------
  # Extract search distance if terra returned it
  # ----------------------------------------------------------
  
  nearest_search_distance <- rep(
    NA_real_,
    nrow(nearest_extract)
  )
  
  
  if (!is.na(distance_column)) {
    
    nearest_search_distance <- suppressWarnings(
      as.numeric(
        nearest_extract[[distance_column]]
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Update the lookup table
  # ----------------------------------------------------------
  
  site_lookup$Vs30_raw[
    rows_to_update
  ] <- nearest_values[
    nearest_source_rows
  ]
  
  
  site_lookup$vs30_raster_cell[
    rows_to_update
  ] <- nearest_cells[
    nearest_source_rows
  ]
  
  
  site_lookup$vs30_cell_x[
    rows_to_update
  ] <- nearest_x[
    nearest_source_rows
  ]
  
  
  site_lookup$vs30_cell_y[
    rows_to_update
  ] <- nearest_y[
    nearest_source_rows
  ]
  
  
  site_lookup$vs30_extraction_source[
    rows_to_update
  ] <- "nearest_nonmissing_raster_cell"
  
  
  # Store the distance returned directly by terra where
  # available. The later distance-calculation block can fill
  # any remaining missing values from the cell coordinates.
  
  if (
    !"vs30_source_distance_m" %in%
    names(site_lookup)
  ) {
    
    site_lookup$vs30_source_distance_m <-
      NA_real_
    
  }
  
  
  site_lookup$vs30_source_distance_m[
    rows_to_update
  ] <- nearest_search_distance[
    nearest_source_rows
  ]
  
  
  cat(
    "Nearest non-missing values assigned:",
    format(
      length(rows_to_update),
      big.mark = ","
    ),
    "\n"
  )
  
  
  cat(
    "Nearest extraction cell column:",
    ifelse(
      is.na(cell_column),
      "not returned; recovered where possible",
      cell_column
    ),
    "\n"
  )
  
}

# ============================================================
# CALCULATE DISTANCE TO SOURCE RASTER-CELL CENTRE
# ============================================================

if (
  !"vs30_source_distance_m" %in%
  names(site_lookup)
) {
  
  site_lookup$vs30_source_distance_m <-
    NA_real_
  
}


valid_source_cell_rows <- which(
  
  !is.na(site_lookup$vs30_cell_x) &
    !is.na(site_lookup$vs30_cell_y) &
    !is.na(site_lookup$Vs30_raw)
  
)


if (length(valid_source_cell_rows) > 0) {
  
  source_cell_points <- terra::vect(
    
    data.frame(
      
      x = site_lookup$vs30_cell_x[
        valid_source_cell_rows
      ],
      
      y = site_lookup$vs30_cell_y[
        valid_source_cell_rows
      ]
      
    ),
    
    geom = c(
      "x",
      "y"
    ),
    
    crs = vs30_crs
    
  )
  
  
  source_distances <- terra::distance(
    
    site_points_raster_crs[
      valid_source_cell_rows
    ],
    
    source_cell_points,
    
    pairwise = TRUE,
    unit = "m"
    
  )
  
  
  site_lookup$vs30_source_distance_m[
    valid_source_cell_rows
  ] <- as.numeric(
    source_distances
  )
  
}


# ============================================================
# CLEAN VS30 VALUES AND CREATE SITE TERMS
# ============================================================

site_lookup <- site_lookup |>
  
  mutate(
    
    vs30_plausible =
      !is.na(Vs30_raw) &
      is.finite(Vs30_raw) &
      Vs30_raw >= minimum_plausible_vs30_m_s &
      Vs30_raw <= maximum_plausible_vs30_m_s,
    
    Vs30_m_s = case_when(
      
      vs30_plausible ~
        Vs30_raw,
      
      TRUE ~
        NA_real_
      
    ),
    
    has_vs30 =
      !is.na(Vs30_m_s),
    
    vs30_nearest_fill_used =
      vs30_extraction_source ==
      "nearest_nonmissing_raster_cell",
    
    ln_Vs30 = case_when(
      
      has_vs30 ~
        log(Vs30_m_s),
      
      TRUE ~
        NA_real_
      
    ),
    
    ln_Vs30_over_760 = case_when(
      
      has_vs30 ~
        log(
          Vs30_m_s /
            vs30_reference_m_s
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    # Vs30-only analytical proxy classes.
    #
    # These should not be interpreted as full regulatory site
    # classifications where other geotechnical criteria apply.
    site_class_vs30_proxy = case_when(
      
      is.na(Vs30_m_s) ~
        NA_character_,
      
      Vs30_m_s > 1500 ~
        "A",
      
      Vs30_m_s >= 760 ~
        "B",
      
      Vs30_m_s >= 360 ~
        "C",
      
      Vs30_m_s >= 180 ~
        "D",
      
      Vs30_m_s < 180 ~
        "E",
      
      TRUE ~
        NA_character_
      
    ),
    
    site_class_broad = case_when(
      
      site_class_vs30_proxy %in%
        c("A", "B") ~
        "rock",
      
      site_class_vs30_proxy == "C" ~
        "stiff_soil",
      
      site_class_vs30_proxy == "D" ~
        "soft_soil",
      
      site_class_vs30_proxy == "E" ~
        "very_soft_soil",
      
      TRUE ~
        NA_character_
      
    ),
    
    vs30_status = case_when(
      
      has_vs30 &
        !vs30_nearest_fill_used ~
        "available_containing_cell",
      
      has_vs30 &
        vs30_nearest_fill_used ~
        "available_nearest_cell",
      
      !is.na(Vs30_raw) &
        !vs30_plausible ~
        "invalid_or_implausible_value",
      
      TRUE ~
        "unavailable"
      
    )
    
  )


# ============================================================
# SAVE UNIQUE SITE LOOKUP
# ============================================================

arrow::write_parquet(
  site_lookup,
  site_lookup_file
)


# ============================================================
# JOIN VS30 VALUES TO ALL FELT RAPID OBSERVATIONS
# ============================================================

site_fields_for_join <- site_lookup |>
  
  select(
    
    report_longitude,
    report_latitude,
    
    site_lookup_id,
    
    Vs30_raw,
    Vs30_m_s,
    
    has_vs30,
    vs30_plausible,
    
    vs30_raster_cell,
    vs30_cell_x,
    vs30_cell_y,
    
    vs30_source_distance_m,
    vs30_extraction_source,
    vs30_nearest_fill_used,
    vs30_status,
    
    ln_Vs30,
    ln_Vs30_over_760,
    
    site_class_vs30_proxy,
    site_class_broad
    
  )


iar_site_all <- iar_rrup |>
  
  left_join(
    
    site_fields_for_join,
    
    by = c(
      "report_longitude",
      "report_latitude"
    )
    
  )


if (nrow(iar_site_all) != nrow(iar_rrup)) {
  
  stop(
    paste0(
      "The Vs30 join changed the number of observations.\n",
      "Before join: ",
      nrow(iar_rrup),
      "\nAfter join: ",
      nrow(iar_site_all)
    ),
    call. = FALSE
  )
  
}


# Observations with invalid coordinates never entered the site
# lookup, so assign explicit status values.

iar_site_all <- iar_site_all |>
  
  mutate(
    
    has_vs30 = replace_na(
      has_vs30,
      FALSE
    ),
    
    vs30_plausible = replace_na(
      vs30_plausible,
      FALSE
    ),
    
    vs30_nearest_fill_used = replace_na(
      vs30_nearest_fill_used,
      FALSE
    ),
    
    vs30_status = case_when(
      
      !valid_site_coordinates ~
        "unavailable_invalid_coordinates",
      
      is.na(vs30_status) ~
        "unavailable",
      
      TRUE ~
        vs30_status
      
    ),
    
    site_model_eligible =
      has_vs30 &
      vs30_plausible &
      valid_site_coordinates
    
  )


# ============================================================
# CREATE FINAL SITE-MODEL DATASET
# ============================================================

iar_final <- iar_site_all |>
  
  filter(
    site_model_eligible
  ) |>
  
  arrange(
    publicid,
    Repi_km,
    cell_id
  )


if (nrow(iar_final) == 0) {
  
  stop(
    paste0(
      "No observations have valid Vs30 values.\n",
      "Check the raster file, raster layer, CRS and ",
      "Felt RAPID coordinates."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# DIAGNOSTICS
# ============================================================

vs30_extraction_summary <- iar_site_all |>
  
  count(
    vs30_status,
    name = "observation_count"
  ) |>
  
  mutate(
    
    percentage =
      100 *
      observation_count /
      sum(observation_count)
    
  ) |>
  
  arrange(
    desc(observation_count)
  )


vs30_unique_site_summary <- site_lookup |>
  
  count(
    vs30_status,
    name = "unique_site_count"
  ) |>
  
  mutate(
    
    percentage =
      100 *
      unique_site_count /
      sum(unique_site_count)
    
  ) |>
  
  arrange(
    desc(unique_site_count)
  )


vs30_site_class_summary <- iar_final |>
  
  count(
    site_class_vs30_proxy,
    site_class_broad,
    name = "observation_count"
  ) |>
  
  mutate(
    
    percentage =
      100 *
      observation_count /
      sum(observation_count)
    
  ) |>
  
  arrange(
    site_class_vs30_proxy
  )


vs30_missing_or_invalid <- iar_site_all |>
  
  filter(
    !site_model_eligible
  ) |>
  
  select(
    
    publicid,
    cell_id = any_of("cell_id"),
    
    report_longitude,
    report_latitude,
    
    reported_mmi,
    report_count = any_of("report_count"),
    
    valid_site_coordinates,
    
    Vs30_raw,
    Vs30_m_s,
    
    vs30_status,
    vs30_extraction_source,
    vs30_source_distance_m
    
  )


raster_extent <- terra::ext(
  vs30_raster
)


vs30_raster_metadata <- tibble(
  
  source_file =
    normalizePath(
      vs30_input_file,
      winslash = "/",
      mustWork = TRUE
    ),
  
  selected_layer_index =
    selected_layer_index,
  
  selected_layer_name =
    available_layer_names[
      selected_layer_index
    ],
  
  raster_rows =
    terra::nrow(vs30_raster),
  
  raster_columns =
    terra::ncol(vs30_raster),
  
  raster_cells =
    terra::ncell(vs30_raster),
  
  raster_resolution_x =
    terra::res(vs30_raster)[[1]],
  
  raster_resolution_y =
    terra::res(vs30_raster)[[2]],
  
  xmin =
    raster_extent$xmin,
  
  xmax =
    raster_extent$xmax,
  
  ymin =
    raster_extent$ymin,
  
  ymax =
    raster_extent$ymax,
  
  crs =
    vs30_crs,
  
  extraction_method =
    vs30_extraction_method,
  
  nearest_search_distance_m =
    maximum_nearest_search_distance_m,
  
  minimum_plausible_vs30_m_s =
    minimum_plausible_vs30_m_s,
  
  maximum_plausible_vs30_m_s =
    maximum_plausible_vs30_m_s
  
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

# Complete Script 6 dataset with Vs30 fields, including rows
# with missing or invalid site values.

arrow::write_parquet(
  iar_site_all,
  all_output_file
)


# Final site-model dataset containing valid Vs30 values.

arrow::write_parquet(
  iar_final,
  final_output_file
)


readr::write_csv(
  vs30_extraction_summary,
  file.path(
    diagnostic_dir,
    "vs30_extraction_summary.csv"
  )
)


readr::write_csv(
  vs30_unique_site_summary,
  file.path(
    diagnostic_dir,
    "vs30_unique_site_summary.csv"
  )
)


readr::write_csv(
  vs30_site_class_summary,
  file.path(
    diagnostic_dir,
    "vs30_site_class_summary.csv"
  )
)


readr::write_csv(
  vs30_missing_or_invalid,
  file.path(
    diagnostic_dir,
    "vs30_missing_or_invalid_observations.csv"
  )
)


readr::write_csv(
  vs30_raster_metadata,
  file.path(
    diagnostic_dir,
    "vs30_raster_metadata.csv"
  )
)


# ============================================================
# VERIFY OUTPUTS
# ============================================================

output_files <- c(
  
  site_lookup_file,
  all_output_file,
  final_output_file,
  
  file.path(
    diagnostic_dir,
    "vs30_extraction_summary.csv"
  ),
  
  file.path(
    diagnostic_dir,
    "vs30_unique_site_summary.csv"
  ),
  
  file.path(
    diagnostic_dir,
    "vs30_site_class_summary.csv"
  ),
  
  file.path(
    diagnostic_dir,
    "vs30_missing_or_invalid_observations.csv"
  ),
  
  file.path(
    diagnostic_dir,
    "vs30_raster_metadata.csv"
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
      "One or more Script 7 outputs were not written correctly:\n",
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
  "\nVs30 site-condition processing complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Input observations:",
  format(
    nrow(iar_rrup),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique site locations:",
  format(
    nrow(site_lookup),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique sites with valid Vs30:",
  format(
    sum(
      site_lookup$has_vs30,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique sites using nearest-cell fill:",
  format(
    sum(
      site_lookup$vs30_nearest_fill_used,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Observations with valid Vs30:",
  format(
    nrow(iar_final),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Observations without valid Vs30:",
  format(
    nrow(iar_site_all) -
      nrow(iar_final),
    big.mark = ","
  ),
  "\n"
)


if (nrow(iar_final) > 0) {
  
  cat(
    "Vs30 range:",
    round(
      min(
        iar_final$Vs30_m_s,
        na.rm = TRUE
      ),
      1
    ),
    "to",
    round(
      max(
        iar_final$Vs30_m_s,
        na.rm = TRUE
      ),
      1
    ),
    "m/s\n"
  )
  
  cat(
    "Median Vs30:",
    round(
      median(
        iar_final$Vs30_m_s,
        na.rm = TRUE
      ),
      1
    ),
    "m/s\n"
  )
  
}


cat(
  "\nUnique-site lookup:\n",
  normalizePath(
    site_lookup_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nComplete dataset with site-condition fields:\n",
  normalizePath(
    all_output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nFinal dataset with valid Vs30:\n",
  normalizePath(
    final_output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)