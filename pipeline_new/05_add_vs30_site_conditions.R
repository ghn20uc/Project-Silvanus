# Extract Vs30 at each Felt RAPID reporting location.
# Save the Foster et al. NZ Vs30 GeoTIFF as data_raw/vs30/combined_vs30.tif.

library(tidyverse)
library(arrow)

vs30_input_file <- file.path("data_raw", "vs30", "combined_vs30.tif")
vs30_layer <- 1L
maximum_nearest_search_distance_m <- 2000
minimum_plausible_vs30_m_s <- 50
maximum_plausible_vs30_m_s <- 3000

input_file <- file.path("data_processed", "NZ_FR_IAR_dataset.parquet")
site_lookup_file <- file.path("data_processed", "vs30_site_lookup.parquet")
all_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_with_site_conditions.parquet"
)
final_output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_final_dataset.parquet"
)

if (!file.exists(input_file)) {
  stop("Modelling dataset not found. Run Script 04 first.")
}

if (!file.exists(vs30_input_file)) {
  stop("Vs30 raster not found at ", vs30_input_file, ".")
}

iar <- read_parquet(input_file) |>
  as_tibble() |>
  mutate(
    valid_site_coordinates =
      !is.na(report_longitude) &
      !is.na(report_latitude) &
      between(report_longitude, -180, 180) &
      between(report_latitude, -90, 90)
  )

vs30_raster <- terra::rast(vs30_input_file)

if (is.character(vs30_layer)) {
  if (!vs30_layer %in% names(vs30_raster)) {
    stop("Vs30 layer not found: ", vs30_layer)
  }
  vs30_raster <- vs30_raster[[vs30_layer]]
} else {
  if (vs30_layer < 1 || vs30_layer > terra::nlyr(vs30_raster)) {
    stop("vs30_layer is outside the available raster layers.")
  }
  vs30_raster <- vs30_raster[[vs30_layer]]
}

names(vs30_raster) <- "Vs30_raw"

if (is.na(terra::crs(vs30_raster)) || terra::crs(vs30_raster) == "") {
  stop("The Vs30 raster does not have a coordinate reference system.")
}

site_lookup <- iar |>
  filter(valid_site_coordinates) |>
  distinct(report_longitude, report_latitude) |>
  arrange(report_longitude, report_latitude) |>
  mutate(site_lookup_id = row_number())

if (nrow(site_lookup) == 0) {
  stop("No valid Felt RAPID coordinates are available for Vs30 extraction.")
}

site_points <- terra::vect(
  site_lookup,
  geom = c("report_longitude", "report_latitude"),
  crs = "EPSG:4326"
) |>
  terra::project(terra::crs(vs30_raster))

extracted <- terra::extract(
  vs30_raster,
  site_points,
  method = "simple",
  cells = TRUE,
  ID = TRUE,
  search_radius = maximum_nearest_search_distance_m
) |>
  as_tibble() |>
  arrange(ID)

if (nrow(extracted) != nrow(site_lookup)) {
  stop("Vs30 extraction did not return one result for each site.")
}

site_lookup <- site_lookup |>
  mutate(
    Vs30_raw = suppressWarnings(as.numeric(extracted$Vs30_raw)),
    vs30_raster_cell = suppressWarnings(as.numeric(extracted$cell)),
    vs30_plausible =
      !is.na(Vs30_raw) &
      between(
        Vs30_raw,
        minimum_plausible_vs30_m_s,
        maximum_plausible_vs30_m_s
      ),
    Vs30_m_s = if_else(vs30_plausible, Vs30_raw, NA_real_),
    has_vs30 = !is.na(Vs30_m_s)
  )

write_parquet(site_lookup, site_lookup_file)

iar_with_vs30 <- iar |>
  left_join(
    site_lookup,
    by = c("report_longitude", "report_latitude")
  ) |>
  mutate(
    has_vs30 = replace_na(has_vs30, FALSE),
    vs30_plausible = replace_na(vs30_plausible, FALSE)
  )

iar_final <- iar_with_vs30 |>
  filter(has_vs30) |>
  arrange(publicid, Repi_km, cell_id)

if (nrow(iar_final) == 0) {
  stop("No modelling observations have a valid Vs30 value.")
}

write_parquet(iar_with_vs30, all_output_file)
write_parquet(iar_final, final_output_file)

message("Saved ", nrow(iar_final), " observations with Vs30 values")
