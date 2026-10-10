# 05_add_cs30_site_conditions
# Thomas Rautao 2026

# Add Vs30 site-condition values to the Felt RAPID analysis dataset.

# Median Vs30 raster Foster et al. (2019)
# https://onlinelibrary.wiley.com/doi/10.1193/121118EQS281M

# Supplementary file esp4bf02915-sup-0003.tif must be downloaded from the above
# source and saved as data_raw/vs30/foster_2019_vs30_median.tif. 

# Vs30 is extracted once for each unique Felt-cell coordinate and then joined
# back to every earthquake-cell observation at that location.

library(tidyverse)
library(arrow)

vs30_input_file <- file.path(
  "data_raw",
  "vs30",
  "foster_2019_vs30_median.tif"
)

# When the raster cell directly beneath a reporting point is NA, look for
# for the nearest non-missing raster cell within this distance. The extracted
# value remains NA when no suitable cell is found.
maximum_nearest_search_distance_m <- 2000

# Remove obvious vs30 errors
minimum_plausible_vs30_m_s <- 50
maximum_plausible_vs30_m_s <- 3000

input_file <- file.path("data_processed", "NZ_FR_IAR_dataset.parquet")
output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_final_dataset.parquet"
)

if (!file.exists(input_file)) {
  stop("Analysis dataset not found. Run Script 04 first.")
}

if (!file.exists(vs30_input_file)) {
  stop("Vs30 raster not found at ", vs30_input_file, ".")
}

iar <- read_parquet(input_file) |>
  as_tibble()

required_columns <- c("report_longitude", "report_latitude")
missing_columns <- setdiff(required_columns, names(iar))

if (length(missing_columns) > 0) {
  stop("Missing site columns: ", paste(missing_columns, collapse = ", "))
}

# The 2019 median-Vs30 GeoTIFF contains one raster layer
vs30_raster <- terra::rast(vs30_input_file)

if (terra::nlyr(vs30_raster) != 1) {
  stop("The Foster et al. (2019) median-Vs30 file should contain one layer.")
}

names(vs30_raster) <- "Vs30_raw"

# The paper raster uses New Zealand Transverse Mercator 2000. EPSG:2193.
terra::crs(vs30_raster) <- "EPSG:2193"

site_lookup <- iar |>
  filter(
    !is.na(report_longitude),
    !is.na(report_latitude),
    between(report_longitude, -180, 180),
    between(report_latitude, -90, 90)
  ) |>
  distinct(report_longitude, report_latitude) |>
  arrange(report_longitude, report_latitude) |>
  mutate(site_id = row_number())

if (nrow(site_lookup) == 0) {
  stop("No valid Felt RAPID coordinates are available for Vs30 extraction.")
}

# Convert the raster extent to longitude/latitude before projecting the Felt
# sites. Some Felt reports can be outside the NZ bounds and sending there 
# coordinates to the NZTM transformation can produce a "Point outside of projection
# domain" warning.
vs30_extent_wgs84 <- terra::project(
  terra::ext(vs30_raster),
  terra::crs(vs30_raster),
  "EPSG:4326"
)

sites_for_extraction <- site_lookup |>
  filter(
    between(
      report_longitude,
      terra::xmin(vs30_extent_wgs84),
      terra::xmax(vs30_extent_wgs84)
    ),
    between(
      report_latitude,
      terra::ymin(vs30_extent_wgs84),
      terra::ymax(vs30_extent_wgs84)
    )
  )

if (nrow(sites_for_extraction) == 0) {
  stop("No Felt RAPID coordinates overlap the Vs30 raster extent.")
}

# The remaining Felt coordinates are within the geographic area covered by the
# raster. Projecting the points rather than the raster preserves the original
# Vs30 grid and values and also makes the search radius use metres.
site_points <- terra::vect(
  sites_for_extraction,
  geom = c("report_longitude", "report_latitude"),
  crs = "EPSG:4326"
) |>
  terra::project(terra::crs(vs30_raster))

projected_site_ids <- terra::values(site_points)$site_id

extracted <- terra::extract(
  vs30_raster,
  site_points,
  method = "simple",
  cells = TRUE,
  ID = TRUE,
  search_radius = maximum_nearest_search_distance_m
) |>
  as_tibble()

if (
  any(is.na(extracted$ID)) ||
  any(!extracted$ID %in% seq_along(projected_site_ids)) ||
  anyDuplicated(extracted$ID) > 0
) {
  stop("Vs30 extraction returned invalid or duplicated point identifiers.")
}

# Values outside the plausible range become NA.
extracted_values <- extracted |>
  transmute(
    site_id = projected_site_ids[ID],
    Vs30_raw = suppressWarnings(as.numeric(Vs30_raw))
  )

site_lookup <- site_lookup |>
  left_join(extracted_values, by = "site_id") |>
  mutate(
    Vs30_m_s = if_else(
      !is.na(Vs30_raw) &
        between(
          Vs30_raw,
          minimum_plausible_vs30_m_s,
          maximum_plausible_vs30_m_s
        ),
      Vs30_raw,
      NA_real_
    )
  ) |>
  select(report_longitude, report_latitude, Vs30_m_s)

iar_final <- iar |>
  left_join(
    site_lookup,
    by = c("report_longitude", "report_latitude")
  ) |>
  arrange(publicid, Repi_km, cell_id)

# Rows without Vs30 are retained with Vs30_m_s = NA
valid_vs30_count <- sum(!is.na(iar_final$Vs30_m_s))

if (valid_vs30_count == 0) {
  stop("Vs30 extraction did not produce any usable values.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_parquet(iar_final, output_file)

message(
  "Saved ", nrow(iar_final), " observations; ",
  valid_vs30_count, " have a valid Vs30 value"
)
