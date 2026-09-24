# 07_build_clean_modelling_dataset.R
# DATA309 - create an improved, auditable modelling dataset
#
# IMPORTANT DESIGN CHOICES:
# 1) Keep the original Script 06 output unchanged.
# 2) Use Moment Tensor Mw only for the strict modelling dataset.
# 3) Do NOT yet delete earthquakes because they have <=5 MMI cells or few Felt reports.
#    Instead, add quality flags and sensitivity tables so the threshold can be justified.
# 4) Do NOT change earthquakes with depth == 33 km. Flag them for investigation.
# 5) Remove the Vs30/760 reference transform from the new analysis datasets.
#    Raw Vs30 and ln(Vs30) are retained.
# 6) Keep the API reported_mmi as the response for now, while adding reconstructed
#    mean/median/mode/skewness fields for auditing and the central-tendency discussion.

library(tidyverse)
library(arrow)

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------

site_path <- "data_processed/NZ_FR_IAR_with_site_conditions.parquet"
counts_path <- "data_raw/felt_rapid_mmi_counts.parquet"
catalogue_path <- "data_processed/catalogue_with_felt_mt_strong.parquet"

out_dir <- "data_processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

quality_path <- file.path(out_dir, "NZ_FR_IAR_quality_flagged.parquet")
strict_path <- file.path(out_dir, "NZ_FR_IAR_modelling_strict_mtMw.parquet")
quality_summary_path <- file.path(out_dir, "DATA309_data_quality_summary.csv")
filter_sensitivity_path <- file.path(out_dir, "DATA309_filter_sensitivity.csv")
event_quality_path <- file.path(out_dir, "DATA309_event_quality.csv")

# -------------------------------------------------------------------
# 1. Read pipeline outputs
# -------------------------------------------------------------------

site <- read_parquet(site_path) |> as_tibble()
mmi_counts <- read_parquet(counts_path) |> as_tibble()
catalogue <- read_parquet(catalogue_path) |> as_tibble()

cat("\n=== INPUT ===\n")
cat("Cell observations:", nrow(site), "\n")
cat("Earthquakes in cell dataset:", n_distinct(site$publicid), "\n")

# -------------------------------------------------------------------
# 2. Helper functions for MMI central tendency
# -------------------------------------------------------------------

weighted_median_discrete <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- as.integer(round(w[ok]))

  if (length(x) == 0 || sum(w) == 0) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  cw <- cumsum(w)
  n <- sum(w)

  value_at_position <- function(pos) {
    x[which(cw >= pos)[1]]
  }

  if (n %% 2 == 1) {
    value_at_position((n + 1) / 2)
  } else {
    mean(c(
      value_at_position(n / 2),
      value_at_position(n / 2 + 1)
    ))
  }
}

weighted_skewness <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]

  if (length(x) < 2 || sum(w) <= 0) {
    return(NA_real_)
  }

  mu <- weighted.mean(x, w)
  m2 <- sum(w * (x - mu)^2) / sum(w)

  if (!is.finite(m2) || m2 <= 0) {
    return(NA_real_)
  }

  m3 <- sum(w * (x - mu)^3) / sum(w)
  m3 / (m2^(3 / 2))
}

# -------------------------------------------------------------------
# 3. Reconstruct cell-level mean / median / mode / skewness
#    These are audit variables; reported_mmi remains the response for now.
# -------------------------------------------------------------------

cell_central <- mmi_counts |>
  filter(
    is.finite(mmi_level),
    is.finite(report_count_at_mmi),
    report_count_at_mmi > 0
  ) |>
  group_by(publicid, cell_id) |>
  summarise(
    reports_from_mmi_counts = sum(report_count_at_mmi),
    reconstructed_mean_mmi =
      weighted.mean(mmi_level, report_count_at_mmi),
    reconstructed_median_mmi =
      weighted_median_discrete(mmi_level, report_count_at_mmi),
    max_mode_count = max(report_count_at_mmi),
    n_modes = sum(report_count_at_mmi == max_mode_count),
    reconstructed_mode_mmi = if_else(
      n_modes == 1L,
      mmi_level[which.max(report_count_at_mmi)],
      NA_real_
    ),
    reconstructed_skewness =
      weighted_skewness(mmi_level, report_count_at_mmi),
    .groups = "drop"
  ) |>
  mutate(
    mode_is_tied = n_modes > 1
  )

# -------------------------------------------------------------------
# 4. Bring in strict Moment Tensor Mw and catalogue audit fields
# -------------------------------------------------------------------

catalogue_key <- catalogue |>
  transmute(
    publicid,
    mt_Mw_from_catalogue_table = mt_Mw,
    catalogue_depth_km = depth,
    catalogue_felt_report_count = felt_report_count
  ) |>
  distinct(publicid, .keep_all = TRUE)

quality <- site |>
  left_join(catalogue_key, by = "publicid") |>
  left_join(cell_central, by = c("publicid", "cell_id"))

# If Script 06 already contains mt_Mw, prefer it; otherwise use the joined field.
if ("mt_Mw" %in% names(quality)) {
  quality <- quality |>
    mutate(Mw_strict = coalesce(mt_Mw, mt_Mw_from_catalogue_table))
} else {
  quality <- quality |>
    mutate(Mw_strict = mt_Mw_from_catalogue_table)
}

# Event depth used by the attenuation dataset.
if ("event_depth_km" %in% names(quality)) {
  quality <- quality |>
    mutate(depth_for_model_km = event_depth_km)
} else {
  quality <- quality |>
    mutate(depth_for_model_km = catalogue_depth_km)
}

# -------------------------------------------------------------------
# 5. Event-level quality information
# -------------------------------------------------------------------

event_quality <- quality |>
  group_by(publicid) |>
  summarise(
    event_mmi_cells = n_distinct(cell_id),
    event_reports_from_cells = sum(report_count, na.rm = TRUE),
    catalogue_felt_report_count = first(catalogue_felt_report_count),
    Mw_strict = first(Mw_strict),
    depth_for_model_km = first(depth_for_model_km),
    .groups = "drop"
  ) |>
  mutate(
    flag_event_4_or_fewer_cells = event_mmi_cells <= 4,
    flag_event_5_or_fewer_cells = event_mmi_cells <= 5,
    flag_event_5_or_fewer_reports = event_reports_from_cells <= 5,
    flag_event_10_or_fewer_reports = event_reports_from_cells <= 10,
    flag_depth_exactly_33km =
      !is.na(depth_for_model_km) &
      abs(depth_for_model_km - 33) < 1e-10,
    flag_missing_mt_Mw = is.na(Mw_strict)
  )

write_csv(event_quality, event_quality_path)

quality <- quality |>
  left_join(
    event_quality |>
      select(
        publicid,
        event_mmi_cells,
        event_reports_from_cells,
        flag_event_4_or_fewer_cells,
        flag_event_5_or_fewer_cells,
        flag_event_5_or_fewer_reports,
        flag_event_10_or_fewer_reports,
        flag_depth_exactly_33km,
        flag_missing_mt_Mw
      ),
    by = "publicid"
  )

# -------------------------------------------------------------------
# 6. Cell-level quality flags
# -------------------------------------------------------------------

quality <- quality |>
  mutate(
    # Report-count sensitivity flags: do not filter yet.
    flag_cell_only_1_report = report_count == 1,
    flag_cell_fewer_than_2_reports = report_count < 2,
    flag_cell_fewer_than_3_reports = report_count < 3,
    flag_cell_fewer_than_5_reports = report_count < 5,

    # Vs30 provenance.
    flag_vs30_nearest_cell =
      vs30_status == "available_nearest_cell",
    flag_vs30_unavailable =
      is.na(Vs30_m_s) |
      vs30_status == "unavailable",

    # Central-tendency audit.
    api_vs_reconstructed_median_difference =
      reported_mmi - reconstructed_median_mmi,
    flag_api_median_mismatch =
      !is.na(reconstructed_median_mmi) &
      abs(api_vs_reconstructed_median_difference) > 1e-10,

    # Core modelling completeness.
    flag_core_model_missing =
      is.na(reported_mmi) |
      is.na(Mw_strict) |
      is.na(Rhypo_km) |
      Rhypo_km <= 0 |
      is.na(Vs30_m_s) |
      Vs30_m_s <= 0 |
      is.na(depth_for_model_km)
  )

# Raw Vs30 plus an unreferenced log transform.
quality <- quality |>
  mutate(
    ln_Vs30_clean = if_else(
      !is.na(Vs30_m_s) & Vs30_m_s > 0,
      log(Vs30_m_s),
      NA_real_
    )
  )

# -------------------------------------------------------------------
# 7. Remove fields Hazel has said should not be used in the final analysis
# -------------------------------------------------------------------

# Remove weighted-mean fields if present.
weighted_mean_cols <- names(quality)[str_detect(
  names(quality),
  regex("weighted.*mean|mean.*weighted", ignore_case = TRUE)
)]

# Remove Vs30 reference-value fields if present.
vs30_reference_cols <- intersect(
  c(
    "ln_Vs30_over_760",
    "vs30_reference_m_s",
    "Vs30_reference_m_s"
  ),
  names(quality)
)

# Remove the temporary joined MT field; Mw_strict is the analysis field.
temporary_cols <- intersect(
  c("mt_Mw_from_catalogue_table"),
  names(quality)
)

quality <- quality |>
  select(
    -any_of(weighted_mean_cols),
    -any_of(vs30_reference_cols),
    -any_of(temporary_cols)
  )

# -------------------------------------------------------------------
# 8. Save the NON-DESTRUCTIVE quality-flagged dataset
# -------------------------------------------------------------------

write_parquet(quality, quality_path)

# -------------------------------------------------------------------
# 9. Create a strict base modelling dataset
#
# This implements Hazel's "Moment Tensor Mw" instruction and only removes
# observations that cannot be used for the core Mw + Rhypo + Vs30 model.
#
# It deliberately does NOT:
#   - require >=5 cells per earthquake
#   - require >=5 reports per cell/event
#   - require Mw >= 5
#   - remove depth == 33 km
#
# Those decisions remain sensitivity/meeting decisions.
# -------------------------------------------------------------------

strict_model <- quality |>
  filter(!flag_core_model_missing) |>
  mutate(
    Mw = Mw_strict
  )

write_parquet(strict_model, strict_path)

# -------------------------------------------------------------------
# 10. Sensitivity table: show what candidate thresholds would cost
# -------------------------------------------------------------------

base_n_rows <- nrow(strict_model)
base_n_events <- n_distinct(strict_model$publicid)

cell_report_thresholds <- c(1, 2, 3, 5, 10)
event_cell_thresholds <- c(1, 3, 5, 10, 20)
event_report_thresholds <- c(1, 3, 5, 10, 20, 50)

cell_sensitivity <- map_dfr(cell_report_thresholds, function(t) {
  d <- strict_model |>
    filter(report_count >= t)

  tibble(
    filter_type = "minimum reports per cell",
    threshold = t,
    rows_retained = nrow(d),
    pct_rows_retained = 100 * nrow(d) / base_n_rows,
    earthquakes_retained = n_distinct(d$publicid),
    pct_earthquakes_retained =
      100 * n_distinct(d$publicid) / base_n_events
  )
})

event_cell_sensitivity <- map_dfr(event_cell_thresholds, function(t) {
  d <- strict_model |>
    filter(event_mmi_cells >= t)

  tibble(
    filter_type = "minimum MMI cells per earthquake",
    threshold = t,
    rows_retained = nrow(d),
    pct_rows_retained = 100 * nrow(d) / base_n_rows,
    earthquakes_retained = n_distinct(d$publicid),
    pct_earthquakes_retained =
      100 * n_distinct(d$publicid) / base_n_events
  )
})

event_report_sensitivity <- map_dfr(event_report_thresholds, function(t) {
  d <- strict_model |>
    filter(event_reports_from_cells >= t)

  tibble(
    filter_type = "minimum Felt reports per earthquake",
    threshold = t,
    rows_retained = nrow(d),
    pct_rows_retained = 100 * nrow(d) / base_n_rows,
    earthquakes_retained = n_distinct(d$publicid),
    pct_earthquakes_retained =
      100 * n_distinct(d$publicid) / base_n_events
  )
})

filter_sensitivity <- bind_rows(
  cell_sensitivity,
  event_cell_sensitivity,
  event_report_sensitivity
)

write_csv(filter_sensitivity, filter_sensitivity_path)

# -------------------------------------------------------------------
# 11. Quality summary
# -------------------------------------------------------------------

quality_summary <- tibble(
  metric = c(
    "original_cell_rows",
    "original_earthquakes",
    "strict_model_cell_rows",
    "strict_model_earthquakes",
    "strict_model_row_retention_pct",
    "rows_missing_mt_Mw",
    "rows_missing_Vs30",
    "rows_using_nearest_Vs30",
    "events_depth_exactly_33km",
    "events_5_or_fewer_MMI_cells",
    "cells_with_api_reconstructed_median_mismatch",
    "cells_with_tied_mode"
  ),
  value = c(
    nrow(quality),
    n_distinct(quality$publicid),
    nrow(strict_model),
    n_distinct(strict_model$publicid),
    100 * nrow(strict_model) / nrow(quality),
    sum(is.na(quality$Mw_strict)),
    sum(is.na(quality$Vs30_m_s)),
    sum(quality$flag_vs30_nearest_cell, na.rm = TRUE),
    sum(event_quality$flag_depth_exactly_33km, na.rm = TRUE),
    sum(event_quality$flag_event_5_or_fewer_cells, na.rm = TRUE),
    sum(quality$flag_api_median_mismatch, na.rm = TRUE),
    sum(quality$mode_is_tied, na.rm = TRUE)
  )
)

write_csv(quality_summary, quality_summary_path)

# -------------------------------------------------------------------
# 12. Console report
# -------------------------------------------------------------------

cat("\n=== IMPROVED DATASET COMPLETE ===\n")
cat("Original rows:", nrow(quality), "\n")
cat("Original earthquakes:", n_distinct(quality$publicid), "\n")
cat("Strict model rows:", nrow(strict_model), "\n")
cat("Strict model earthquakes:", n_distinct(strict_model$publicid), "\n")
cat(
  "Strict row retention:",
  round(100 * nrow(strict_model) / nrow(quality), 2),
  "%\n"
)

cat("\nFiles created:\n")
cat("  ", quality_path, "\n")
cat("  ", strict_path, "\n")
cat("  ", event_quality_path, "\n")
cat("  ", quality_summary_path, "\n")
cat("  ", filter_sensitivity_path, "\n")

cat("\nIMPORTANT:\n")
cat("- No event/cell-count threshold has been forced yet.\n")
cat("- No minimum MT Mw threshold has been forced yet.\n")
cat("- Depth == 33 km has only been flagged, not altered.\n")
cat("- The strict dataset uses Moment Tensor Mw only.\n")
