# Add GeoNet centroid moment-tensor data to the Felt RAPID catalogue.

library(tidyverse)
library(arrow)

input_file <- file.path(
  "data_processed",
  "catalogue_with_felt_summary.parquet"
)
moment_tensor_file <- file.path("data_raw", "moment_tensor.parquet")
output_file <- file.path(
  "data_processed",
  "catalogue_with_felt_mt.parquet"
)

mt_url <- paste0(
  "https://raw.githubusercontent.com/GeoNet/data/main/",
  "moment-tensor/GeoNet_CMT_solutions.csv"
)

if (!file.exists(input_file)) {
  stop("Felt RAPID catalogue not found. Run Script 02 first.")
}

catalogue <- read_parquet(input_file) |>
  as_tibble()

mt_raw <- read_csv(
  mt_url,
  col_types = cols(.default = col_character()),
  na = c("", "NA", "N/A", "NaN", "NULL", "null", "-"),
  show_col_types = FALSE,
  progress = FALSE
) |>
  rename_with(~ str_to_lower(str_trim(.x)))

if ("event_id" %in% names(mt_raw) && !"publicid" %in% names(mt_raw)) {
  mt_raw <- rename(mt_raw, publicid = event_id)
}

required_columns <- c(
  "publicid", "date", "mw", "cd", "ns", "vr",
  "strike1", "dip1", "rake1", "strike2", "dip2", "rake2"
)
missing_columns <- setdiff(required_columns, names(mt_raw))

if (length(missing_columns) > 0) {
  stop("Missing moment-tensor columns: ", paste(missing_columns, collapse = ", "))
}

numeric_columns <- intersect(
  c(
    "latitude", "longitude", "strike1", "dip1", "rake1",
    "strike2", "dip2", "rake2", "ml", "mw", "mo", "cd",
    "ns", "dc", "vr", "mxx", "mxy", "mxz", "myy", "myz", "mzz"
  ),
  names(mt_raw)
)

parse_solution_time <- function(x) {
  x <- str_remove(str_trim(as.character(x)), "\\.0+$")
  suppressWarnings(as.POSIXct(x, format = "%Y%m%d%H%M%S", tz = "UTC"))
}

fault_style_from_rake <- function(rake) {
  rake <- ((rake + 180) %% 360) - 180

  case_when(
    is.na(rake) ~ NA_character_,
    abs(rake) <= 30 | abs(abs(rake) - 180) <= 30 ~ "strike-slip",
    rake > 30 & rake < 150 ~ "reverse",
    rake < -30 & rake > -150 ~ "normal",
    TRUE ~ "oblique"
  )
}

mt <- mt_raw |>
  mutate(
    source_row = row_number(),
    publicid = na_if(str_trim(publicid), ""),
    solution_time = parse_solution_time(date),
    across(all_of(numeric_columns), ~ suppressWarnings(parse_double(.x)))
  ) |>
  filter(!is.na(publicid), publicid != "9999999") |>
  transmute(
    publicid,
    source_row,
    solution_time,
    Mw = if_else(between(mw, 0, 10), mw, NA_real_),
    mt_centroid_latitude = latitude,
    mt_centroid_longitude = longitude,
    mt_centroid_depth_km = if_else(between(cd, 0, 700), cd, NA_real_),
    mt_station_count = as.integer(round(ns)),
    mt_variance_reduction_percent = vr,
    mt_double_couple_percent = dc,
    mt_scalar_moment_dyne_cm = mo,
    mt_strike1 = strike1,
    mt_dip1 = dip1,
    mt_rake1 = rake1,
    mt_strike2 = strike2,
    mt_dip2 = dip2,
    mt_rake2 = rake2,
    mt_method = method
  ) |>
  mutate(
    fault_style_np1 = fault_style_from_rake(mt_rake1),
    fault_style_np2 = fault_style_from_rake(mt_rake2),
    fault_style = if_else(
      !is.na(fault_style_np1) & fault_style_np1 == fault_style_np2,
      fault_style_np1,
      NA_character_
    ),
    completeness = rowSums(
      across(
        c(
          Mw, mt_centroid_depth_km, mt_strike1, mt_dip1, mt_rake1,
          mt_strike2, mt_dip2, mt_rake2
        ),
        ~ !is.na(.x)
      )
    )
  )

# Prefer a complete solution with Mw, then use fit, station count and date.
mt_preferred <- mt |>
  arrange(
    publicid,
    desc(!is.na(Mw)),
    desc(completeness),
    desc(coalesce(mt_variance_reduction_percent, -Inf)),
    desc(coalesce(as.numeric(mt_station_count), -Inf)),
    desc(solution_time),
    desc(source_row)
  ) |>
  distinct(publicid, .keep_all = TRUE) |>
  select(-source_row, -completeness) |>
  mutate(
    has_moment_tensor = TRUE,
    Mw_source = if_else(!is.na(Mw), "GeoNet_CMT", NA_character_)
  )

if (nrow(mt_preferred) == 0) {
  stop("No moment-tensor solutions had a usable GeoNet public ID.")
}

dir.create(dirname(moment_tensor_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

write_parquet(mt_preferred, moment_tensor_file)

catalogue_with_mt <- catalogue |>
  select(-any_of(c("Mw", "Mw_model", "Mw_source", "fault_style"))) |>
  left_join(mt_preferred, by = "publicid") |>
  mutate(has_moment_tensor = replace_na(has_moment_tensor, FALSE))

write_parquet(catalogue_with_mt, output_file)

message(
  "Matched moment tensors to ",
  sum(catalogue_with_mt$has_moment_tensor),
  " of ", nrow(catalogue_with_mt),
  " catalogue earthquakes"
)
