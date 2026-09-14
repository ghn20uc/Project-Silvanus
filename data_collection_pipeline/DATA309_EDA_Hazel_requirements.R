# DATA309 EDA
# Run this script from the Data_Collection_Pipeline project directory.
# Expected input files:
#   data_processed/NZ_FR_IAR_with_site_conditions.parquet
#   data_raw/felt_rapid_mmi_counts.parquet
#   data_raw/felt_rapid_observations.parquet
#   data_processed/catalogue_with_felt_mt_strong.parquet
#
# Outputs are written to: eda_outputs/

library(tidyverse)
library(arrow)

out_dir <- "eda_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# -------------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------------

site <- read_parquet(
  "data_processed/NZ_FR_IAR_with_site_conditions.parquet"
) |>
  as_tibble()

mmi_counts <- read_parquet(
  "data_raw/felt_rapid_mmi_counts.parquet"
) |>
  as_tibble()

felt_obs <- read_parquet(
  "data_raw/felt_rapid_observations.parquet"
) |>
  as_tibble()

catalogue <- read_parquet(
  "data_processed/catalogue_with_felt_mt_strong.parquet"
) |>
  as_tibble()


cat("\n=== DATA LOADED ===\n")
cat("Site-condition rows:", nrow(site), "\n")
cat("MMI-count rows:", nrow(mmi_counts), "\n")
cat("Felt observation rows:", nrow(felt_obs), "\n")
cat("Catalogue earthquakes:", nrow(catalogue), "\n")


# Events represented in the modelling/site dataset
final_event_ids <- site |>
  distinct(publicid)


events_final <- catalogue |>
  semi_join(
    final_event_ids,
    by = "publicid"
  ) |>
  distinct(
    publicid,
    .keep_all = TRUE
  )


events_with_felt <- catalogue |>
  filter(
    has_felt_rapid %in% TRUE
  ) |>
  distinct(
    publicid,
    .keep_all = TRUE
  )


overview <- tibble(
  metric = c(
    "catalogue_events",
    "events_with_felt_rapid",
    "events_in_site_dataset",
    "felt_cells_in_site_dataset",
    "unique_felt_sites",
    "mmi_count_rows"
  ),
  value = c(
    nrow(catalogue),
    nrow(events_with_felt),
    nrow(events_final),
    nrow(site),
    nrow(
      site |>
        distinct(
          report_longitude,
          report_latitude
        )
    ),
    nrow(mmi_counts)
  )
)


write_csv(
  overview,
  file.path(
    out_dir,
    "00_dataset_overview.csv"
  )
)

print(overview)


# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------

num_summary <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      tibble(
        n = 0,
        missing = NA_integer_,
        min = NA_real_,
        q1 = NA_real_,
        median = NA_real_,
        mean = NA_real_,
        q3 = NA_real_,
        max = NA_real_,
        sd = NA_real_
      )
    )
  }
  
  tibble(
    n = length(x),
    min = min(x),
    q1 = unname(
      quantile(
        x,
        0.25
      )
    ),
    median = median(x),
    mean = mean(x),
    q3 = unname(
      quantile(
        x,
        0.75
      )
    ),
    max = max(x),
    sd = sd(x)
  )
}


weighted_median_discrete <- function(x, w) {
  
  ok <-
    is.finite(x) &
    is.finite(w) &
    w > 0
  
  x <- x[ok]
  w <- as.integer(
    round(
      w[ok]
    )
  )
  
  if (
    length(x) == 0 ||
    sum(w) == 0
  ) {
    return(
      NA_real_
    )
  }
  
  ord <- order(x)
  
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w)
  n <- sum(w)
  
  
  value_at_position <- function(pos) {
    
    x[
      which(
        cw >= pos
      )[1]
    ]
  }
  
  
  if (n %% 2 == 1) {
    
    value_at_position(
      (n + 1) / 2
    )
    
  } else {
    
    mean(
      c(
        value_at_position(
          n / 2
        ),
        value_at_position(
          n / 2 + 1
        )
      )
    )
  }
}


weighted_skewness <- function(x, w) {
  
  ok <-
    is.finite(x) &
    is.finite(w) &
    w > 0
  
  x <- x[ok]
  w <- w[ok]
  
  if (
    length(x) < 2 ||
    sum(w) <= 0
  ) {
    return(
      NA_real_
    )
  }
  
  mu <- weighted.mean(
    x,
    w
  )
  
  m2 <-
    sum(
      w *
        (x - mu)^2
    ) /
    sum(w)
  
  
  if (
    !is.finite(m2) ||
    m2 <= 0
  ) {
    return(
      NA_real_
    )
  }
  
  m3 <-
    sum(
      w *
        (x - mu)^3
    ) /
    sum(w)
  
  m3 /
    (
      m2^(3 / 2)
    )
}


save_plot <- function(
    p,
    filename,
    width = 8,
    height = 5.5
) {
  
  ggsave(
    filename = file.path(
      out_dir,
      filename
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}


# -------------------------------------------------------------------
# 2. Moment Tensor Mw distribution
# IMPORTANT: use mt_Mw, not M_model fallback
# -------------------------------------------------------------------

mw_summary <- events_final |>
  summarise(
    final_events = n(),
    
    events_with_mt_Mw =
      sum(
        !is.na(mt_Mw)
      ),
    
    pct_with_mt_Mw =
      100 *
      mean(
        !is.na(mt_Mw)
      ),
    
    min_mt_Mw =
      min(
        mt_Mw,
        na.rm = TRUE
      ),
    
    median_mt_Mw =
      median(
        mt_Mw,
        na.rm = TRUE
      ),
    
    mean_mt_Mw =
      mean(
        mt_Mw,
        na.rm = TRUE
      ),
    
    max_mt_Mw =
      max(
        mt_Mw,
        na.rm = TRUE
      )
  )


write_csv(
  mw_summary,
  file.path(
    out_dir,
    "01_mt_Mw_summary.csv"
  )
)

print(mw_summary)


p_mw <- events_final |>
  filter(
    is.finite(mt_Mw)
  ) |>
  ggplot(
    aes(
      x = mt_Mw
    )
  ) +
  geom_histogram(
    binwidth = 0.2,
    boundary = 0
  ) +
  labs(
    title =
      "Moment Tensor Mw distribution",
    subtitle =
      "One row per earthquake in the site-condition dataset",
    x =
      "Moment magnitude (mt_Mw)",
    y =
      "Number of earthquakes"
  ) +
  theme_minimal()


save_plot(
  p_mw,
  "01_mt_Mw_distribution.png"
)


# Candidate Mw thresholds:
# RETENTION ONLY, not a recommendation
mw_thresholds <- c(
  5.0,
  5.25,
  5.5,
  5.75,
  6.0,
  6.5,
  7.0
)


mw_threshold_retention <- map_dfr(
  mw_thresholds,
  function(t) {
    
    retained_ids <- events_final |>
      filter(
        is.finite(mt_Mw),
        mt_Mw >= t
      ) |>
      pull(
        publicid
      )
    
    
    tibble(
      minimum_mt_Mw = t,
      
      earthquakes_retained =
        length(
          unique(
            retained_ids
          )
        ),
      
      cells_retained =
        sum(
          site$publicid %in%
            retained_ids
        ),
      
      pct_final_earthquakes_retained =
        100 *
        length(
          unique(
            retained_ids
          )
        ) /
        nrow(events_final),
      
      pct_cells_retained =
        100 *
        sum(
          site$publicid %in%
            retained_ids
        ) /
        nrow(site)
    )
  }
)


write_csv(
  mw_threshold_retention,
  file.path(
    out_dir,
    "01b_mt_Mw_threshold_retention.csv"
  )
)


# -------------------------------------------------------------------
# 3. Earthquake depth distribution
# -------------------------------------------------------------------

depth_summary <- events_final |>
  summarise(
    n_events = n(),
    
    missing_depth =
      sum(
        is.na(depth)
      ),
    
    min_depth_km =
      min(
        depth,
        na.rm = TRUE
      ),
    
    q1_depth_km =
      unname(
        quantile(
          depth,
          0.25,
          na.rm = TRUE
        )
      ),
    
    median_depth_km =
      median(
        depth,
        na.rm = TRUE
      ),
    
    mean_depth_km =
      mean(
        depth,
        na.rm = TRUE
      ),
    
    q3_depth_km =
      unname(
        quantile(
          depth,
          0.75,
          na.rm = TRUE
        )
      ),
    
    max_depth_km =
      max(
        depth,
        na.rm = TRUE
      )
  )


write_csv(
  depth_summary,
  file.path(
    out_dir,
    "02_depth_summary.csv"
  )
)

print(depth_summary)


p_depth <- events_final |>
  filter(
    is.finite(depth)
  ) |>
  ggplot(
    aes(
      x = depth
    )
  ) +
  geom_histogram(
    bins = 35
  ) +
  labs(
    title =
      "Earthquake depth distribution",
    subtitle =
      "One row per earthquake",
    x =
      "Depth (km)",
    y =
      "Number of earthquakes"
  ) +
  theme_minimal()


save_plot(
  p_depth,
  "02_depth_distribution.png"
)


# -------------------------------------------------------------------
# 4. Felt report-count distributions
# Cell level and earthquake level are intentionally separated
# -------------------------------------------------------------------

cell_report_summary <- felt_obs |>
  summarise(
    cells = n(),
    
    missing_report_count =
      sum(
        is.na(report_count)
      ),
    
    min =
      min(
        report_count,
        na.rm = TRUE
      ),
    
    q1 =
      unname(
        quantile(
          report_count,
          0.25,
          na.rm = TRUE
        )
      ),
    
    median =
      median(
        report_count,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        report_count,
        na.rm = TRUE
      ),
    
    q3 =
      unname(
        quantile(
          report_count,
          0.75,
          na.rm = TRUE
        )
      ),
    
    max =
      max(
        report_count,
        na.rm = TRUE
      )
  )


write_csv(
  cell_report_summary,
  file.path(
    out_dir,
    "03_report_count_per_cell_summary.csv"
  )
)

print(cell_report_summary)


p_cell_reports <- felt_obs |>
  filter(
    is.finite(report_count),
    report_count > 0
  ) |>
  ggplot(
    aes(
      x = report_count
    )
  ) +
  geom_histogram(
    bins = 40
  ) +
  scale_x_log10() +
  labs(
    title =
      "Felt report count per reporting cell",
    x =
      "Reports per cell (log10 scale)",
    y =
      "Number of cells"
  ) +
  theme_minimal()


save_plot(
  p_cell_reports,
  "03_report_count_per_cell_distribution.png"
)


event_report_summary <- events_final |>
  summarise(
    earthquakes = n(),
    
    missing_felt_report_count =
      sum(
        is.na(
          felt_report_count
        )
      ),
    
    min =
      min(
        felt_report_count,
        na.rm = TRUE
      ),
    
    q1 =
      unname(
        quantile(
          felt_report_count,
          0.25,
          na.rm = TRUE
        )
      ),
    
    median =
      median(
        felt_report_count,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        felt_report_count,
        na.rm = TRUE
      ),
    
    q3 =
      unname(
        quantile(
          felt_report_count,
          0.75,
          na.rm = TRUE
        )
      ),
    
    max =
      max(
        felt_report_count,
        na.rm = TRUE
      )
  )


write_csv(
  event_report_summary,
  file.path(
    out_dir,
    "04_felt_report_count_per_earthquake_summary.csv"
  )
)

print(event_report_summary)


p_event_reports <- events_final |>
  filter(
    is.finite(
      felt_report_count
    ),
    felt_report_count > 0
  ) |>
  ggplot(
    aes(
      x = felt_report_count
    )
  ) +
  geom_histogram(
    bins = 40
  ) +
  scale_x_log10() +
  labs(
    title =
      "Total Felt reports per earthquake",
    x =
      "Felt reports per earthquake (log10 scale)",
    y =
      "Number of earthquakes"
  ) +
  theme_minimal()


save_plot(
  p_event_reports,
  "04_felt_report_count_per_earthquake.png"
)


# Report-count threshold retention at CELL level
cell_thresholds <- c(
  1,
  2,
  3,
  5,
  6,
  10,
  20
)


cell_threshold_retention <- map_dfr(
  cell_thresholds,
  function(t) {
    
    d <- felt_obs |>
      filter(
        is.finite(
          report_count
        ),
        report_count >= t
      )
    
    
    tibble(
      minimum_reports_per_cell = t,
      
      cells_retained =
        nrow(d),
      
      earthquakes_retained =
        n_distinct(
          d$publicid
        ),
      
      pct_cells_retained =
        100 *
        nrow(d) /
        nrow(felt_obs)
    )
  }
)


write_csv(
  cell_threshold_retention,
  file.path(
    out_dir,
    "04b_cell_report_threshold_retention.csv"
  )
)


# -------------------------------------------------------------------
# 5. MMI cells per earthquake
# -------------------------------------------------------------------

cells_per_event <- site |>
  distinct(
    publicid,
    cell_id
  ) |>
  count(
    publicid,
    name = "mmi_cells"
  )


cells_per_event_summary <- cells_per_event |>
  summarise(
    earthquakes = n(),
    
    min =
      min(
        mmi_cells
      ),
    
    q1 =
      unname(
        quantile(
          mmi_cells,
          0.25
        )
      ),
    
    median =
      median(
        mmi_cells
      ),
    
    mean =
      mean(
        mmi_cells
      ),
    
    q3 =
      unname(
        quantile(
          mmi_cells,
          0.75
        )
      ),
    
    max =
      max(
        mmi_cells
      ),
    
    events_with_4_or_fewer_cells =
      sum(
        mmi_cells <= 4
      ),
    
    pct_with_4_or_fewer_cells =
      100 *
      mean(
        mmi_cells <= 4
      ),
    
    events_with_5_or_fewer_cells =
      sum(
        mmi_cells <= 5
      ),
    
    pct_with_5_or_fewer_cells =
      100 *
      mean(
        mmi_cells <= 5
      )
  )


write_csv(
  cells_per_event_summary,
  file.path(
    out_dir,
    "05_mmi_cells_per_earthquake_summary.csv"
  )
)


write_csv(
  cells_per_event,
  file.path(
    out_dir,
    "05b_mmi_cells_per_earthquake.csv"
  )
)


print(
  cells_per_event_summary
)


p_cells_event <- cells_per_event |>
  ggplot(
    aes(
      x = mmi_cells
    )
  ) +
  geom_histogram(
    bins = 40
  ) +
  scale_x_log10() +
  labs(
    title =
      "MMI reporting cells per earthquake",
    x =
      "Number of cells (log10 scale)",
    y =
      "Number of earthquakes"
  ) +
  theme_minimal()


save_plot(
  p_cells_event,
  "05_mmi_cells_per_earthquake.png"
)


# -------------------------------------------------------------------
# 6. Mean / median / mode and skewness for each reporting cell
# Uses the MMI-category count table, not the pipeline's weighted-mean field.
# -------------------------------------------------------------------

cell_central <- mmi_counts |>
  filter(
    is.finite(
      mmi_level
    ),
    is.finite(
      report_count_at_mmi
    ),
    report_count_at_mmi > 0
  ) |>
  group_by(
    publicid,
    cell_id
  ) |>
  summarise(
    total_reports_from_counts =
      sum(
        report_count_at_mmi
      ),
    
    mean_mmi =
      weighted.mean(
        mmi_level,
        report_count_at_mmi
      ),
    
    median_mmi =
      weighted_median_discrete(
        mmi_level,
        report_count_at_mmi
      ),
    
    max_mode_count =
      max(
        report_count_at_mmi
      ),
    
    n_modes =
      sum(
        report_count_at_mmi ==
          max_mode_count
      ),
    
    mode_mmi =
      if_else(
        n_modes == 1L,
        mmi_level[
          which.max(
            report_count_at_mmi
          )
        ],
        NA_real_
      ),
    
    skewness =
      weighted_skewness(
        mmi_level,
        report_count_at_mmi
      ),
    
    .groups = "drop"
  ) |>
  left_join(
    felt_obs |>
      select(
        publicid,
        cell_id,
        reported_mmi,
        report_count
      ) |>
      distinct(),
    by = c(
      "publicid",
      "cell_id"
    )
  ) |>
  mutate(
    mean_minus_median =
      mean_mmi -
      median_mmi,
    
    api_median_minus_reconstructed =
      reported_mmi -
      median_mmi,
    
    mode_is_tied =
      n_modes > 1
  )


write_csv(
  cell_central,
  file.path(
    out_dir,
    "06_cell_central_tendency_and_skewness.csv"
  )
)


central_summary <- cell_central |>
  summarise(
    cells = n(),
    
    cells_with_unique_mode =
      sum(
        !is.na(
          mode_mmi
        )
      ),
    
    pct_tied_mode =
      100 *
      mean(
        mode_is_tied
      ),
    
    mean_abs_mean_median_difference =
      mean(
        abs(
          mean_minus_median
        ),
        na.rm = TRUE
      ),
    
    pct_mean_equals_median =
      100 *
      mean(
        abs(
          mean_minus_median
        ) < 1e-12,
        na.rm = TRUE
      ),
    
    pct_api_median_matches_reconstructed =
      100 *
      mean(
        abs(
          api_median_minus_reconstructed
        ) < 1e-12,
        na.rm = TRUE
      )
  )


write_csv(
  central_summary,
  file.path(
    out_dir,
    "06b_central_tendency_summary.csv"
  )
)


print(
  central_summary
)


p_mean_median <- cell_central |>
  filter(
    is.finite(
      mean_mmi
    ),
    is.finite(
      median_mmi
    )
  ) |>
  ggplot(
    aes(
      x = median_mmi,
      y = mean_mmi
    )
  ) +
  geom_jitter(
    width = 0.08,
    height = 0.08,
    alpha = 0.2
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  labs(
    title =
      "Cell-level MMI: mean versus median",
    x =
      "Median MMI",
    y =
      "Mean MMI"
  ) +
  theme_minimal()


save_plot(
  p_mean_median,
  "06_mean_vs_median_MMI.png"
)


p_mode_median <- cell_central |>
  filter(
    is.finite(
      mode_mmi
    ),
    is.finite(
      median_mmi
    )
  ) |>
  ggplot(
    aes(
      x = median_mmi,
      y = mode_mmi
    )
  ) +
  geom_jitter(
    width = 0.08,
    height = 0.08,
    alpha = 0.2
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  labs(
    title =
      "Cell-level MMI: mode versus median",
    subtitle =
      "Cells with tied modes are excluded",
    x =
      "Median MMI",
    y =
      "Unique mode MMI"
  ) +
  theme_minimal()


save_plot(
  p_mode_median,
  "06_mode_vs_median_MMI.png"
)


# -------------------------------------------------------------------
# 6b. Report-count threshold and within-cell MMI stability
#
# Purpose:
# Assess whether requiring more Felt reports per cell improves the
# stability/precision of the cell-level MMI estimate, and quantify
# the amount of modelling/site data lost at each candidate threshold.
#
# Note:
# within_cell_mmi_sd measures disagreement among reports in a cell.
# estimated_mean_se is used only as a precision proxy and should not
# be interpreted as the uncertainty of the GeoNet reported_mmi itself.
#
# Only cells represented in the site-condition / modelling dataset
# are retained for this analysis.
# -------------------------------------------------------------------

cell_stability <- mmi_counts |>
  filter(
    is.finite(
      mmi_level
    ),
    is.finite(
      report_count_at_mmi
    ),
    report_count_at_mmi > 0
  ) |>
  group_by(
    publicid,
    cell_id
  ) |>
  summarise(
    total_reports =
      sum(
        report_count_at_mmi
      ),
    
    weighted_mean_mmi =
      weighted.mean(
        mmi_level,
        report_count_at_mmi
      ),
    
    within_cell_mmi_sd = {
      
      n_total <-
        sum(
          report_count_at_mmi
        )
      
      mu <-
        weighted.mean(
          mmi_level,
          report_count_at_mmi
        )
      
      
      if (n_total > 1) {
        
        sqrt(
          sum(
            report_count_at_mmi *
              (
                mmi_level -
                  mu
              )^2
          ) /
            (
              n_total - 1
            )
        )
        
      } else {
        
        NA_real_
        
      }
    },
    
    .groups = "drop"
  ) |>
  mutate(
    estimated_mean_se =
      within_cell_mmi_sd /
      sqrt(
        total_reports
      )
  ) |>
  inner_join(
    site |>
      select(
        publicid,
        cell_id,
        reported_mmi,
        report_count
      ) |>
      distinct(),
    by = c(
      "publicid",
      "cell_id"
    )
  )


# Candidate minimum report-count thresholds
stability_thresholds <- c(
  1,
  2,
  3,
  5,
  6,
  10,
  20
)


threshold_stability <- map_dfr(
  stability_thresholds,
  function(t) {
    
    d <- cell_stability |>
      filter(
        is.finite(
          report_count
        ),
        report_count >= t
      )
    
    
    tibble(
      minimum_reports_per_cell = t,
      
      cells_retained =
        nrow(d),
      
      pct_cells_retained =
        100 *
        nrow(d) /
        nrow(
          cell_stability
        ),
      
      earthquakes_retained =
        n_distinct(
          d$publicid
        ),
      
      median_within_cell_sd =
        median(
          d$within_cell_mmi_sd,
          na.rm = TRUE
        ),
      
      mean_within_cell_sd =
        mean(
          d$within_cell_mmi_sd,
          na.rm = TRUE
        ),
      
      median_estimated_mean_se =
        median(
          d$estimated_mean_se,
          na.rm = TRUE
        ),
      
      mean_estimated_mean_se =
        mean(
          d$estimated_mean_se,
          na.rm = TRUE
        )
    )
  }
)


write_csv(
  threshold_stability,
  file.path(
    out_dir,
    "06c_report_threshold_stability.csv"
  )
)


print(
  threshold_stability
)


# Data retention plot

p_threshold_retention <- threshold_stability |>
  ggplot(
    aes(
      x =
        minimum_reports_per_cell,
      y =
        pct_cells_retained
    )
  ) +
  geom_line() +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks =
      stability_thresholds
  ) +
  labs(
    title =
      "Data retained under minimum report-count thresholds",
    x =
      "Minimum Felt reports per cell",
    y =
      "Cells retained (%)"
  ) +
  theme_minimal()


save_plot(
  p_threshold_retention,
  "06c_report_threshold_retention.png"
)


# Precision proxy plot

p_threshold_precision <- threshold_stability |>
  ggplot(
    aes(
      x =
        minimum_reports_per_cell,
      y =
        median_estimated_mean_se
    )
  ) +
  geom_line() +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks =
      stability_thresholds
  ) +
  labs(
    title =
      "Estimated cell-level MMI precision by report threshold",
    subtitle =
      "Lower median SE indicates greater sampling precision",
    x =
      "Minimum Felt reports per cell",
    y =
      "Median estimated SE of within-cell mean MMI"
  ) +
  theme_minimal()


save_plot(
  p_threshold_precision,
  "06d_report_threshold_precision.png"
)


# -------------------------------------------------------------------
# 7. Does skewness change as median MMI increases?
# -------------------------------------------------------------------

skew_by_median <- cell_central |>
  filter(
    is.finite(
      median_mmi
    ),
    is.finite(
      skewness
    )
  ) |>
  group_by(
    median_mmi
  ) |>
  summarise(
    cells = n(),
    
    mean_skewness =
      mean(
        skewness
      ),
    
    median_skewness =
      median(
        skewness
      ),
    
    q1_skewness =
      unname(
        quantile(
          skewness,
          0.25
        )
      ),
    
    q3_skewness =
      unname(
        quantile(
          skewness,
          0.75
        )
      ),
    
    .groups = "drop"
  )


write_csv(
  skew_by_median,
  file.path(
    out_dir,
    "07_skewness_by_median_MMI.csv"
  )
)


p_skew <- cell_central |>
  filter(
    is.finite(
      median_mmi
    ),
    is.finite(
      skewness
    )
  ) |>
  ggplot(
    aes(
      x =
        factor(
          median_mmi
        ),
      y =
        skewness
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.15
  ) +
  labs(
    title =
      "Skewness of Felt-report distributions by median MMI",
    x =
      "Median MMI",
    y =
      "Weighted skewness"
  ) +
  theme_minimal()


save_plot(
  p_skew,
  "07_skewness_vs_median_MMI.png"
)


# -------------------------------------------------------------------
# 8. Core distributions: reported MMI, Rhypo, Vs30
# -------------------------------------------------------------------

p_mmi <- site |>
  filter(
    is.finite(
      reported_mmi
    )
  ) |>
  ggplot(
    aes(
      x =
        reported_mmi
    )
  ) +
  geom_bar() +
  labs(
    title =
      "Reported MMI distribution",
    x =
      "Reported MMI",
    y =
      "Number of reporting cells"
  ) +
  theme_minimal()


save_plot(
  p_mmi,
  "08_reported_MMI_distribution.png"
)


p_rhypo <- site |>
  filter(
    is.finite(
      Rhypo_km
    ),
    Rhypo_km > 0
  ) |>
  ggplot(
    aes(
      x =
        Rhypo_km
    )
  ) +
  geom_histogram(
    bins = 40
  ) +
  scale_x_log10() +
  labs(
    title =
      "Hypocentral distance distribution",
    x =
      "Rhypo (km, log10 scale)",
    y =
      "Number of reporting cells"
  ) +
  theme_minimal()


save_plot(
  p_rhypo,
  "09_Rhypo_distribution.png"
)


p_vs30 <- site |>
  filter(
    is.finite(
      Vs30_m_s
    ),
    Vs30_m_s > 0
  ) |>
  ggplot(
    aes(
      x =
        Vs30_m_s
    )
  ) +
  geom_histogram(
    bins = 40
  ) +
  labs(
    title =
      "Vs30 distribution",
    x =
      "Vs30 (m/s)",
    y =
      "Number of reporting cells"
  ) +
  theme_minimal()


save_plot(
  p_vs30,
  "10_Vs30_distribution.png"
)


# Vs30 extraction status

if (
  "vs30_status" %in%
  names(site)
) {
  
  vs30_status_summary <- site |>
    count(
      vs30_status,
      sort = TRUE
    ) |>
    mutate(
      pct =
        100 *
        n /
        sum(n)
    )
  
  
  write_csv(
    vs30_status_summary,
    file.path(
      out_dir,
      "10b_vs30_status_summary.csv"
    )
  )
}


# -------------------------------------------------------------------
# 9. Core relationships for attenuation EDA
# -------------------------------------------------------------------

# Ensure strict moment-tensor Mw is available in the cell-level data

if (
  !"mt_Mw" %in%
  names(site)
) {
  
  site <- site |>
    left_join(
      catalogue |>
        select(
          publicid,
          mt_Mw
        ),
      by = "publicid"
    )
}


p_mmi_distance <- site |>
  filter(
    is.finite(
      reported_mmi
    ),
    is.finite(
      Rhypo_km
    ),
    Rhypo_km > 0
  ) |>
  mutate(
    log10_Rhypo =
      log10(
        Rhypo_km
      )
  ) |>
  ggplot(
    aes(
      x =
        log10_Rhypo,
      y =
        reported_mmi
    )
  ) +
  geom_jitter(
    height = 0.08,
    width = 0,
    alpha = 0.12
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    title =
      "Reported MMI versus hypocentral distance",
    x =
      "log10(Rhypo km)",
    y =
      "Reported MMI"
  ) +
  theme_minimal()


save_plot(
  p_mmi_distance,
  "11_MMI_vs_Rhypo.png"
)


p_mmi_mw <- site |>
  filter(
    is.finite(
      reported_mmi
    ),
    is.finite(
      mt_Mw
    )
  ) |>
  ggplot(
    aes(
      x =
        mt_Mw,
      y =
        reported_mmi
    )
  ) +
  geom_jitter(
    width = 0.02,
    height = 0.08,
    alpha = 0.12
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    title =
      "Reported MMI versus Moment Tensor Mw",
    x =
      "Moment Tensor Mw",
    y =
      "Reported MMI"
  ) +
  theme_minimal()


save_plot(
  p_mmi_mw,
  "12_MMI_vs_mt_Mw.png"
)


p_mmi_vs30 <- site |>
  filter(
    is.finite(
      reported_mmi
    ),
    is.finite(
      Vs30_m_s
    ),
    Vs30_m_s > 0
  ) |>
  ggplot(
    aes(
      x =
        log(
          Vs30_m_s
        ),
      y =
        reported_mmi
    )
  ) +
  geom_jitter(
    height = 0.08,
    width = 0,
    alpha = 0.12
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    title =
      "Reported MMI versus Vs30",
    x =
      "ln(Vs30)",
    y =
      "Reported MMI"
  ) +
  theme_minimal()


save_plot(
  p_mmi_vs30,
  "13_MMI_vs_Vs30.png"
)


# -------------------------------------------------------------------
# 10. Missingness of key modelling fields
# -------------------------------------------------------------------

key_fields <- intersect(
  c(
    "reported_mmi",
    "mt_Mw",
    "Rhypo_km",
    "event_depth_km",
    "Vs30_m_s",
    "report_count"
  ),
  names(site)
)


missing_summary <- map_dfr(
  key_fields,
  function(v) {
    
    tibble(
      variable = v,
      
      missing_n =
        sum(
          is.na(
            site[[v]]
          )
        ),
      
      missing_pct =
        100 *
        mean(
          is.na(
            site[[v]]
          )
        )
    )
  }
)


write_csv(
  missing_summary,
  file.path(
    out_dir,
    "14_key_variable_missingness.csv"
  )
)


# -------------------------------------------------------------------
# 11. Short console report
# -------------------------------------------------------------------

cat(
  "\n\n==============================\n"
)

cat(
  "EDA COMPLETE\n"
)

cat(
  "==============================\n"
)

cat(
  "Outputs written to:",
  normalizePath(
    out_dir
  ),
  "\n\n"
)


cat(
  "Key files to inspect first:\n"
)

cat(
  "  01_mt_Mw_summary.csv\n"
)

cat(
  "  02_depth_summary.csv\n"
)

cat(
  "  04b_cell_report_threshold_retention.csv\n"
)

cat(
  "  05_mmi_cells_per_earthquake_summary.csv\n"
)

cat(
  "  06b_central_tendency_summary.csv\n"
)

cat(
  "  06c_report_threshold_stability.csv\n"
)

cat(
  "  07_skewness_by_median_MMI.csv\n"
)

cat(
  "  10b_vs30_status_summary.csv (if created)\n\n"
)


cat(
  "Key figures:\n"
)

cat(
  "  01_mt_Mw_distribution.png\n"
)

cat(
  "  02_depth_distribution.png\n"
)

cat(
  "  03_report_count_per_cell_distribution.png\n"
)

cat(
  "  04_felt_report_count_per_earthquake.png\n"
)

cat(
  "  05_mmi_cells_per_earthquake.png\n"
)

cat(
  "  06_mean_vs_median_MMI.png\n"
)

cat(
  "  06_mode_vs_median_MMI.png\n"
)

cat(
  "  06c_report_threshold_retention.png\n"
)

cat(
  "  06d_report_threshold_precision.png\n"
)

cat(
  "  07_skewness_vs_median_MMI.png\n"
)

cat(
  "  11_MMI_vs_Rhypo.png\n"
)

cat(
  "  12_MMI_vs_mt_Mw.png\n"
)

cat(
  "  13_MMI_vs_Vs30.png\n"
)