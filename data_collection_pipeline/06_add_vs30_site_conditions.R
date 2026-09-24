# Extracts Vs30 site-condition values
# NOTES TO SELF:
# the Foster et al. New Zealand Vs30 GeoTIFF file must be downloaded from 
# https://ucdigitalsms.atlassian.net/wiki/spaces/QuakeCore/pages/3291711635/Vs30+with+QGIS
# and placed in a user-created folder named "vs30" within the data_raw folder
# this is the best source I've found for the "If you want to USE this work, you can download the 
# publication and electronic supplements, including the model itself as GeoTIF format raster files"
# mentioned in the GitHub (https://github.com/fostergeotech/Vs30_NZ), the link to the website
# from the GitHub is defunct. As far as I can tell this is the same file referenced, however we also have the 
# option to make one using the FosterGeotech GitHub code
# Need to decide on what vs30_reference_m_s should be. 760 from https://onlinelibrary.wiley.com/doi/10.1193/063013EQS181M?
# Earthquakes can contain the same Felt RAPID reporting location. Extract Vs30 only once per unique coordinate pair

library(tidyverse)

# Controls

vs30_input_file <- file.path(
  "data_raw",
  "vs30",
  "combined_vs30.tif"
)

# Raster layer to extract.

vs30_layer <- 1L

# Extraction method:

vs30_extraction_method <- "simple"

# Where a point falls in a raster cell with an NA value, search nearby

maximum_nearest_search_distance_m <- 2000

# Plausibility limits

minimum_plausible_vs30_m_s <- 50
maximum_plausible_vs30_m_s <- 3000

# The main continuous site term will be centred on this value

vs30_reference_m_s <- 760

# Set file paths

input_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset.parquet"
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

dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

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

# Load data

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Script 05 output was not found:\n",
      normalizePath(
        input_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 05 first."
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
    "The Script 05 dataset contains zero rows.",
    call. = FALSE
  )
}

# Standardise coordiante types

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

# Load Vs30 Raster

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

# Select raster layer

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

# Validate raster CRS

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

# Build unique site lookup

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

# Create spatial site points

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

# Exact raster cell extraction

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
# Nearest non-missing cell for exact NA values

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

  # Identify returned metadata columns
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

  # Extract returned Vs30 values
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

  # Extract cell numbers when available
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

  # Extract cell-centre coordinates when available
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

  # Recover cell numbers from returned x/y where needed
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

  # Recover cell-centre x/y from cell numbers where needed
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

  # Extract search distance if terra returned it
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

  # Update lookup table
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
}

# Calculate distance to source raster-cell centre

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

# Clean Vs30 values and create site terms

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
# Save unique site lookup

arrow::write_parquet(
  site_lookup,
  site_lookup_file
)

# join Vs30 values to Felt RAPID observations

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

# Assign explicit status values to observations with invalid coordinates

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
# Create final site-model dataset

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

# Save outputs

# Complete Script 06 dataset

arrow::write_parquet(
  iar_site_all,
  all_output_file
)

# Final site-model dataset containing valid Vs30 values

arrow::write_parquet(
  iar_final,
  final_output_file
)

# Completion Summary

cat(
  "\nVs30 site-condition processing complete\n"
)
