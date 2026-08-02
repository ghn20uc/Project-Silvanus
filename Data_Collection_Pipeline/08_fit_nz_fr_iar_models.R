# ============================================================
# 08_fit_nz_fr_iar_models.R
#
# Fit and compare linear mixed-effects intensity attenuation
# models for New Zealand Felt RAPID observations.
#
# Primary response:
#   reported_mmi
#
# Random effect:
#   earthquake-specific intercept, grouped by publicid
#
# Required input:
#   data_processed/NZ_FR_IAR_final_dataset.parquet
#
# Principal outputs:
#   models/NZ_FR_IAR_final_model.rds
#   data_processed/NZ_FR_IAR_model_sample.parquet
#   data_processed/NZ_FR_IAR_with_model_results.parquet
#   model_results/*.csv
#   model_results/plots/*.png
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(arrow)
library(lme4)
library(lmerTest)
library(performance)


# ============================================================
# MODELLING CONTROLS
# ============================================================

# ------------------------------------------------------------
# Analysis period
# ------------------------------------------------------------

analysis_start_date <- as.POSIXct(
  "2012-01-01 00:00:00",
  tz = "UTC"
)


# Leave as NA for no upper limit.

analysis_end_date <- as.POSIXct(
  NA_character_,
  tz = "UTC"
)


# ------------------------------------------------------------
# Magnitude controls
# ------------------------------------------------------------

minimum_model_magnitude <- 2.5


# FALSE:
#   use M_model, which uses recognised Mw where available and
#   otherwise falls back to the GeoNet catalogue magnitude.
#
# TRUE:
#   retain only rows with recognised Mw_model values.
#
# FALSE maximises the usable dataset but means that the
# magnitude predictor is not uniformly Mw.

require_recognised_mw <- FALSE


# Magnitude reference used for centring.

magnitude_reference <- 5.0


# ------------------------------------------------------------
# Distance controls
# ------------------------------------------------------------

# Recommended primary model:
#
#   "Rhypo_km"
#
# Optional finite-rupture sensitivity model:
#
#   "Rrup_approx_preferred_km"
#
# The approximate Rrup field is only available for a subset of
# earthquakes and depends on estimated rupture geometry.

distance_variable <- "Rhypo_km"


minimum_distance_km <- 0
maximum_distance_km <- 500


# Offset prevents log10(0) if approximate Rrup contains zero.

distance_log_offset_km <- 1


# Reference distance used for centring log distance.

distance_reference_km <- 50


# ------------------------------------------------------------
# Depth controls
# ------------------------------------------------------------

minimum_depth_km <- 0
maximum_depth_km <- 700

depth_reference_km <- 15


# ------------------------------------------------------------
# Felt-report controls
# ------------------------------------------------------------

minimum_reports_per_cell <- 1L

minimum_cells_per_event <- 3L

minimum_total_events <- 20L


# Expected MMI range.

minimum_valid_mmi <- 1
maximum_valid_mmi <- 12


# ------------------------------------------------------------
# Weight controls
# ------------------------------------------------------------

# Available methods:
#
#   "none"
#   "sqrt_report_count"
#   "report_count"
#   "log_report_count"
#
# sqrt_report_count is a moderate weighting method that gives
# more influence to cells supported by more reports without
# allowing high-count cells to dominate as strongly as raw
# report-count weights.

weight_method <- "sqrt_report_count"


# Report counts above this value receive the same maximum
# weight. This protects against a small number of very
# high-count cells dominating the model.

maximum_report_count_for_weight <- 25L


# ------------------------------------------------------------
# Optional predictors
# ------------------------------------------------------------

request_tectonic_region_model <- TRUE

request_fault_style_model <- TRUE


# A categorical level must satisfy both thresholds.

minimum_observations_per_factor_level <- 100L

minimum_events_per_factor_level <- 10L


# At least this proportion of the otherwise eligible analysis
# rows must remain after requiring optional factor values.

minimum_optional_complete_fraction <- 0.60


# ------------------------------------------------------------
# Model-selection controls
# ------------------------------------------------------------

# Select the simplest model within this many AIC units of the
# lowest-AIC model.

aic_parsimony_threshold <- 2


# ============================================================
# VALIDATE CONTROLS
# ============================================================

valid_weight_methods <- c(
  "none",
  "sqrt_report_count",
  "report_count",
  "log_report_count"
)


if (!weight_method %in% valid_weight_methods) {
  
  stop(
    paste0(
      "weight_method must be one of:\n",
      paste(
        valid_weight_methods,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


if (
  minimum_distance_km < 0 ||
  maximum_distance_km <= minimum_distance_km
) {
  
  stop(
    "The distance limits are invalid.",
    call. = FALSE
  )
  
}


if (
  minimum_depth_km < 0 ||
  maximum_depth_km <= minimum_depth_km
) {
  
  stop(
    "The depth limits are invalid.",
    call. = FALSE
  )
  
}


if (
  minimum_optional_complete_fraction <= 0 ||
  minimum_optional_complete_fraction > 1
) {
  
  stop(
    paste0(
      "minimum_optional_complete_fraction must be greater ",
      "than zero and no greater than one."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# PATHS
# ============================================================

input_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_final_dataset.parquet"
)

model_directory <- "models"

result_directory <- "model_results"

plot_directory <- file.path(
  result_directory,
  "plots"
)

diagnostic_directory <- file.path(
  "data_processed",
  "diagnostics"
)

model_sample_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_model_sample.parquet"
)

model_results_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_with_model_results.parquet"
)

final_model_file <- file.path(
  model_directory,
  "NZ_FR_IAR_final_model.rds"
)

final_model_ml_file <- file.path(
  model_directory,
  "NZ_FR_IAR_selected_model_ML.rds"
)

candidate_models_file <- file.path(
  model_directory,
  "NZ_FR_IAR_candidate_models_ML.rds"
)


dir.create(
  model_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostic_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# LOAD SCRIPT 7 DATASET
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Script 7 output was not found:\n",
      normalizePath(
        input_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 7 first."
    ),
    call. = FALSE
  )
  
}


iar_source <- arrow::read_parquet(
  input_file
) |>
  as_tibble() |>
  mutate(
    observation_id = row_number()
  )


if (nrow(iar_source) == 0) {
  
  stop(
    "The Script 7 dataset contains zero observations.",
    call. = FALSE
  )
  
}


cat(
  "Script 7 observations loaded:",
  format(
    nrow(iar_source),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "publicid",
  "origintime",
  "reported_mmi",
  "report_count",
  "event_depth_km",
  "Vs30_m_s"
)


missing_required_columns <- setdiff(
  required_columns,
  names(iar_source)
)


if (length(missing_required_columns) > 0) {
  
  stop(
    paste0(
      "The final dataset is missing required columns:\n",
      paste(
        missing_required_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


if (!distance_variable %in% names(iar_source)) {
  
  stop(
    paste0(
      "The selected distance variable was not found:\n",
      distance_variable,
      "\n\nAvailable distance-like columns include:\n",
      paste(
        grep(
          "^(R|distance)",
          names(iar_source),
          value = TRUE,
          ignore.case = TRUE
        ),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# ------------------------------------------------------------
# Add optional columns where absent
# ------------------------------------------------------------

if (!"M_model" %in% names(iar_source)) {
  iar_source$M_model <- NA_real_
}

if (!"Mw_model" %in% names(iar_source)) {
  iar_source$Mw_model <- NA_real_
}

if (!"M_model_is_Mw" %in% names(iar_source)) {
  iar_source$M_model_is_Mw <- !is.na(
    iar_source$Mw_model
  )
}

if (!"tectonic_region" %in% names(iar_source)) {
  iar_source$tectonic_region <- NA_character_
}

if (!"fault_style" %in% names(iar_source)) {
  iar_source$fault_style <- NA_character_
}

if (!"ln_Vs30_over_760" %in% names(iar_source)) {
  
  iar_source$ln_Vs30_over_760 <- case_when(
    
    !is.na(iar_source$Vs30_m_s) &
      iar_source$Vs30_m_s > 0 ~
      log(
        iar_source$Vs30_m_s /
          760
      ),
    
    TRUE ~
      NA_real_
    
  )
  
}


# ============================================================
# TYPE AND LABEL HELPERS
# ============================================================

parse_origin_time <- function(x) {
  
  if (inherits(x, "POSIXt")) {
    
    return(
      as.POSIXct(
        x,
        tz = "UTC"
      )
    )
    
  }
  
  
  suppressWarnings(
    lubridate::ymd_hms(
      x,
      quiet = TRUE,
      tz = "UTC"
    )
  )
  
}


clean_optional_category <- function(x) {
  
  result <- x |>
    as.character() |>
    stringr::str_trim()
  
  
  missing_labels <- c(
    "",
    "NA",
    "N/A",
    "na",
    "n/a",
    "unknown",
    "Unknown",
    "UNKNOWN",
    "unclassified",
    "Unclassified",
    "not_available",
    "not available",
    "none"
  )
  
  
  result[
    result %in% missing_labels
  ] <- NA_character_
  
  
  result
  
}


# ============================================================
# PREPARE MODEL VARIABLES
# ============================================================

iar_prepared <- iar_source |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    origintime = parse_origin_time(
      origintime
    ),
    
    reported_mmi = suppressWarnings(
      as.numeric(reported_mmi)
    ),
    
    report_count = suppressWarnings(
      as.integer(report_count)
    ),
    
    event_depth_km = suppressWarnings(
      as.numeric(event_depth_km)
    ),
    
    Vs30_m_s = suppressWarnings(
      as.numeric(Vs30_m_s)
    ),
    
    M_model = suppressWarnings(
      as.numeric(M_model)
    ),
    
    Mw_model = suppressWarnings(
      as.numeric(Mw_model)
    ),
    
    model_distance_km = suppressWarnings(
      as.numeric(
        .data[[distance_variable]]
      )
    ),
    
    model_magnitude = case_when(
      
      require_recognised_mw ~
        Mw_model,
      
      TRUE ~
        coalesce(
          M_model,
          Mw_model
        )
      
    ),
    
    model_magnitude_is_Mw = case_when(
      
      require_recognised_mw ~
        !is.na(Mw_model),
      
      TRUE ~
        M_model_is_Mw %in% TRUE |
        (
          is.na(M_model) &
            !is.na(Mw_model)
        )
      
    ),
    
    site_term =
      suppressWarnings(
        as.numeric(
          ln_Vs30_over_760
        )
      ),
    
    tectonic_region_model =
      clean_optional_category(
        tectonic_region
      ),
    
    fault_style_model =
      clean_optional_category(
        fault_style
      )
    
  )


# ============================================================
# INITIAL ELIGIBILITY FILTERS
# ============================================================

iar_prepared <- iar_prepared |>
  
  mutate(
    
    initial_exclusion_reason = case_when(
      
      is.na(publicid) ~
        "missing_publicid",
      
      is.na(origintime) ~
        "missing_origintime",
      
      origintime < analysis_start_date ~
        "before_analysis_start_date",
      
      !is.na(analysis_end_date) &
        origintime > analysis_end_date ~
        "after_analysis_end_date",
      
      is.na(reported_mmi) ~
        "missing_mmi",
      
      reported_mmi < minimum_valid_mmi |
        reported_mmi > maximum_valid_mmi ~
        "mmi_outside_valid_range",
      
      is.na(model_magnitude) ~
        "missing_model_magnitude",
      
      model_magnitude <
        minimum_model_magnitude ~
        "below_minimum_model_magnitude",
      
      require_recognised_mw &
        !model_magnitude_is_Mw ~
        "magnitude_not_recognised_Mw",
      
      is.na(model_distance_km) ~
        "missing_model_distance",
      
      model_distance_km <
        minimum_distance_km ~
        "below_minimum_distance",
      
      model_distance_km >
        maximum_distance_km ~
        "beyond_maximum_distance",
      
      is.na(event_depth_km) ~
        "missing_event_depth",
      
      event_depth_km <
        minimum_depth_km |
        event_depth_km >
        maximum_depth_km ~
        "event_depth_outside_valid_range",
      
      is.na(Vs30_m_s) |
        Vs30_m_s <= 0 ~
        "missing_or_invalid_Vs30",
      
      is.na(site_term) ~
        "missing_site_term",
      
      is.na(report_count) ~
        "missing_report_count",
      
      report_count <
        minimum_reports_per_cell ~
        "insufficient_reports_in_cell",
      
      TRUE ~
        NA_character_
      
    )
    
  )


# ------------------------------------------------------------
# Require enough cells within each event
# ------------------------------------------------------------

initial_event_counts <- iar_prepared |>
  
  filter(
    is.na(initial_exclusion_reason)
  ) |>
  
  count(
    publicid,
    name = "eligible_cell_count"
  )


iar_prepared <- iar_prepared |>
  
  left_join(
    initial_event_counts,
    by = "publicid"
  ) |>
  
  mutate(
    
    exclusion_reason = case_when(
      
      !is.na(initial_exclusion_reason) ~
        initial_exclusion_reason,
      
      is.na(eligible_cell_count) |
        eligible_cell_count <
        minimum_cells_per_event ~
        "insufficient_cells_for_event",
      
      TRUE ~
        NA_character_
      
    )
    
  )


analysis_data <- iar_prepared |>
  
  filter(
    is.na(exclusion_reason)
  )


if (nrow(analysis_data) == 0) {
  
  stop(
    "No observations passed the basic modelling filters.",
    call. = FALSE
  )
  
}


if (
  n_distinct(
    analysis_data$publicid
  ) < minimum_total_events
) {
  
  stop(
    paste0(
      "Only ",
      n_distinct(
        analysis_data$publicid
      ),
      " earthquakes passed the basic filters. At least ",
      minimum_total_events,
      " are required by the current controls."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# CENTRE CONTINUOUS PREDICTORS
# ============================================================

analysis_data <- analysis_data |>
  
  mutate(
    
    magnitude_c =
      model_magnitude -
      magnitude_reference,
    
    log_distance =
      log10(
        model_distance_km +
          distance_log_offset_km
      ),
    
    log_distance_reference =
      log10(
        distance_reference_km +
          distance_log_offset_km
      ),
    
    log_distance_c =
      log_distance -
      log_distance_reference,
    
    depth_10km_c =
      (
        event_depth_km -
          depth_reference_km
      ) /
      10
    
  )


# ============================================================
# ASSESS OPTIONAL CATEGORICAL PREDICTORS
# ============================================================

assess_optional_factor <- function(
    data,
    column_name
) {
  
  level_summary <- data |>
    
    filter(
      !is.na(
        .data[[column_name]]
      )
    ) |>
    
    group_by(
      factor_level =
        .data[[column_name]]
    ) |>
    
    summarise(
      
      observation_count =
        n(),
      
      event_count =
        n_distinct(publicid),
      
      .groups = "drop"
      
    ) |>
    
    filter(
      
      observation_count >=
        minimum_observations_per_factor_level,
      
      event_count >=
        minimum_events_per_factor_level
      
    ) |>
    
    arrange(
      desc(observation_count)
    )
  
  
  retained_levels <-
    level_summary$factor_level
  
  
  retained_fraction <- mean(
    data[[column_name]] %in%
      retained_levels
  )
  
  
  list(
    
    usable =
      length(retained_levels) >= 2 &&
      retained_fraction >=
      minimum_optional_complete_fraction,
    
    retained_levels =
      retained_levels,
    
    retained_fraction =
      retained_fraction,
    
    summary =
      level_summary
    
  )
  
}


tectonic_assessment <- assess_optional_factor(
  analysis_data,
  "tectonic_region_model"
)


fault_assessment <- assess_optional_factor(
  analysis_data,
  "fault_style_model"
)


use_tectonic_region <-
  request_tectonic_region_model &&
  tectonic_assessment$usable


use_fault_style <-
  request_fault_style_model &&
  fault_assessment$usable


# ------------------------------------------------------------
# Remove sparse category levels
# ------------------------------------------------------------

if (use_tectonic_region) {
  
  tectonic_values <-
    analysis_data$tectonic_region_model
  
  
  tectonic_values[
    !tectonic_values %in%
      tectonic_assessment$retained_levels
  ] <- NA_character_
  
  
  analysis_data$tectonic_region_model <-
    factor(
      tectonic_values
    )
  
}


if (use_fault_style) {
  
  fault_values <-
    analysis_data$fault_style_model
  
  
  fault_values[
    !fault_values %in%
      fault_assessment$retained_levels
  ] <- NA_character_
  
  
  analysis_data$fault_style_model <-
    factor(
      fault_values
    )
  
}


# ============================================================
# BUILD A COMMON MODEL-COMPARISON SAMPLE
#
# Every candidate model must be fitted to exactly the same rows.
# ============================================================

build_comparison_sample <- function(
    data,
    include_tectonic,
    include_fault
) {
  
  required_model_columns <- c(
    "reported_mmi",
    "publicid",
    "magnitude_c",
    "log_distance_c",
    "depth_10km_c",
    "site_term",
    "report_count"
  )
  
  
  if (include_tectonic) {
    
    required_model_columns <- c(
      required_model_columns,
      "tectonic_region_model"
    )
    
  }
  
  
  if (include_fault) {
    
    required_model_columns <- c(
      required_model_columns,
      "fault_style_model"
    )
    
  }
  
  
  result <- data |>
    
    tidyr::drop_na(
      all_of(
        required_model_columns
      )
    ) |>
    
    group_by(
      publicid
    ) |>
    
    filter(
      n() >= minimum_cells_per_event
    ) |>
    
    ungroup()
  
  
  if (include_tectonic) {
    
    result$tectonic_region_model <-
      droplevels(
        factor(
          result$tectonic_region_model
        )
      )
    
  }
  
  
  if (include_fault) {
    
    result$fault_style_model <-
      droplevels(
        factor(
          result$fault_style_model
        )
      )
    
  }
  
  
  result
  
}


repeat {
  
  comparison_data <- build_comparison_sample(
    
    data = analysis_data,
    
    include_tectonic =
      use_tectonic_region,
    
    include_fault =
      use_fault_style
    
  )
  
  
  retained_fraction <-
    nrow(comparison_data) /
    nrow(analysis_data)
  
  
  too_few_events <-
    n_distinct(
      comparison_data$publicid
    ) <
    minimum_total_events
  
  
  tectonic_has_too_few_levels <-
    use_tectonic_region &&
    nlevels(
      comparison_data$tectonic_region_model
    ) < 2
  
  
  fault_has_too_few_levels <-
    use_fault_style &&
    nlevels(
      comparison_data$fault_style_model
    ) < 2
  
  
  sample_is_too_small <-
    retained_fraction <
    minimum_optional_complete_fraction ||
    too_few_events
  
  
  if (
    !tectonic_has_too_few_levels &&
    !fault_has_too_few_levels &&
    !sample_is_too_small
  ) {
    
    break
    
  }
  
  
  # Remove the least essential optional predictor first.
  
  if (use_fault_style) {
    
    warning(
      paste0(
        "Fault style was removed from the common comparison ",
        "sample because coverage or factor levels were ",
        "insufficient."
      ),
      call. = FALSE
    )
    
    use_fault_style <- FALSE
    
    next
    
  }
  
  
  if (use_tectonic_region) {
    
    warning(
      paste0(
        "Tectonic region was removed from the common ",
        "comparison sample because coverage or factor levels ",
        "were insufficient."
      ),
      call. = FALSE
    )
    
    use_tectonic_region <- FALSE
    
    next
    
  }
  
  
  stop(
    paste0(
      "The base model sample contains insufficient rows or ",
      "earthquakes under the current controls."
    ),
    call. = FALSE
  )
  
}


# ============================================================
# CREATE MODEL WEIGHTS
# ============================================================

comparison_data <- comparison_data |>
  
  mutate(
    
    capped_report_count =
      pmin(
        pmax(
          report_count,
          1L
        ),
        maximum_report_count_for_weight
      ),
    
    raw_model_weight = case_when(
      
      weight_method ==
        "none" ~
        1,
      
      weight_method ==
        "sqrt_report_count" ~
        sqrt(
          capped_report_count
        ),
      
      weight_method ==
        "report_count" ~
        as.numeric(
          capped_report_count
        ),
      
      weight_method ==
        "log_report_count" ~
        log1p(
          capped_report_count
        ),
      
      TRUE ~
        1
      
    ),
    
    # lmer does not normalise prior weights, so explicitly
    # normalise them to mean one.
    
    model_weight =
      raw_model_weight /
      mean(
        raw_model_weight,
        na.rm = TRUE
      )
    
  )


if (
  any(
    !is.finite(
      comparison_data$model_weight
    )
  ) ||
  any(
    comparison_data$model_weight <= 0
  )
) {
  
  stop(
    "Invalid model weights were created.",
    call. = FALSE
  )
  
}


cat(
  "\nCommon model-comparison sample\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Observations:",
  format(
    nrow(comparison_data),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes:",
  format(
    n_distinct(
      comparison_data$publicid
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Distance variable:",
  distance_variable,
  "\n"
)

cat(
  "Magnitude restricted to recognised Mw:",
  require_recognised_mw,
  "\n"
)

cat(
  "Tectonic region included:",
  use_tectonic_region,
  "\n"
)

cat(
  "Fault style included:",
  use_fault_style,
  "\n"
)

cat(
  "Weight method:",
  weight_method,
  "\n"
)


# ============================================================
# DEFINE CANDIDATE MODEL FORMULAS
# ============================================================

base_terms <- c(
  "magnitude_c",
  "log_distance_c",
  "depth_10km_c"
)


model_terms <- list(
  
  M1_Base =
    base_terms,
  
  M2_Site =
    c(
      base_terms,
      "site_term"
    )
  
)


current_terms <- model_terms$M2_Site


if (use_tectonic_region) {
  
  current_terms <- c(
    current_terms,
    "tectonic_region_model"
  )
  
  
  model_terms$M3_Tectonic <-
    current_terms
  
}


if (use_fault_style) {
  
  current_terms <- c(
    current_terms,
    "fault_style_model"
  )
  
  
  model_terms$M4_FaultStyle <-
    current_terms
  
}


model_terms$M5_MagnitudeDistanceInteraction <- c(
  current_terms,
  "magnitude_c:log_distance_c"
)


make_mixed_formula <- function(terms) {
  
  as.formula(
    
    paste(
      
      "reported_mmi ~",
      
      paste(
        terms,
        collapse = " + "
      ),
      
      "+ (1 | publicid)"
      
    )
    
  )
  
}


model_formulas <- purrr::map(
  model_terms,
  make_mixed_formula
)


# ============================================================
# FIT CANDIDATE MODELS USING MAXIMUM LIKELIHOOD
# ============================================================

model_control <- lme4::lmerControl(
  
  optimizer = "bobyqa",
  
  optCtrl = list(
    maxfun = 200000
  )
  
)


fit_candidate_model <- function(model_formula) {
  
  tryCatch(
    
    lmerTest::lmer(
      
      formula =
        model_formula,
      
      data =
        comparison_data,
      
      weights =
        model_weight,
      
      REML =
        FALSE,
      
      control =
        model_control,
      
      na.action =
        na.fail
      
    ),
    
    error = function(e) {
      e
    }
    
  )
  
}


model_attempts <- purrr::map(
  model_formulas,
  fit_candidate_model
)


failed_model_names <- names(
  model_attempts
)[
  purrr::map_lgl(
    model_attempts,
    inherits,
    what = "error"
  )
]


if (length(failed_model_names) > 0) {
  
  warning(
    paste0(
      "The following candidate model(s) failed:\n",
      paste(
        failed_model_names,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


models_ml <- model_attempts[
  
  !purrr::map_lgl(
    model_attempts,
    inherits,
    what = "error"
  )
  
]


if (length(models_ml) == 0) {
  
  stop(
    "All candidate models failed to fit.",
    call. = FALSE
  )
  
}


# ============================================================
# MODEL-COMPARISON TABLE
# ============================================================

extract_convergence_message <- function(model) {
  
  messages <-
    model@optinfo$conv$lme4$messages
  
  
  if (
    is.null(messages) ||
    length(messages) == 0
  ) {
    
    return("ok")
    
  }
  
  
  paste(
    messages,
    collapse = "; "
  )
  
}


model_comparison <- purrr::imap_dfr(
  
  models_ml,
  
  function(model, model_name) {
    
    tibble(
      
      Model =
        model_name,
      
      model_order =
        match(
          model_name,
          names(model_formulas)
        ),
      
      observations =
        stats::nobs(model),
      
      parameters =
        attr(
          logLik(model),
          "df"
        ),
      
      log_likelihood =
        as.numeric(
          logLik(model)
        ),
      
      AIC =
        AIC(model),
      
      BIC =
        BIC(model),
      
      singular =
        lme4::isSingular(
          model,
          tol = 1e-5
        ),
      
      convergence =
        extract_convergence_message(
          model
        ),
      
      formula =
        paste(
          deparse(
            formula(model)
          ),
          collapse = ""
        )
      
    )
    
  }
  
) |>
  
  mutate(
    
    delta_AIC =
      AIC -
      min(
        AIC,
        na.rm = TRUE
      ),
    
    delta_BIC =
      BIC -
      min(
        BIC,
        na.rm = TRUE
      )
    
  ) |>
  
  arrange(
    model_order
  )


print(
  model_comparison |>
    select(
      Model,
      observations,
      parameters,
      AIC,
      delta_AIC,
      BIC,
      delta_BIC,
      singular,
      convergence
    )
)


# ============================================================
# LIKELIHOOD-RATIO TESTS BETWEEN SUCCESSIVE MODELS
# ============================================================

model_names <- names(
  models_ml
)


if (length(model_names) >= 2) {
  
  likelihood_ratio_tests <- purrr::map2_dfr(
    
    model_names[
      seq_len(
        length(model_names) - 1
      )
    ],
    
    model_names[
      2:length(model_names)
    ],
    
    function(reduced_name, full_name) {
      
      reduced_model <-
        models_ml[[reduced_name]]
      
      full_model <-
        models_ml[[full_name]]
      
      
      reduced_loglik <-
        logLik(reduced_model)
      
      full_loglik <-
        logLik(full_model)
      
      
      chi_square <-
        2 *
        (
          as.numeric(full_loglik) -
            as.numeric(reduced_loglik)
        )
      
      
      df_difference <-
        attr(
          full_loglik,
          "df"
        ) -
        attr(
          reduced_loglik,
          "df"
        )
      
      
      tibble(
        
        reduced_model =
          reduced_name,
        
        full_model =
          full_name,
        
        chi_square =
          chi_square,
        
        df_difference =
          df_difference,
        
        p_value = case_when(
          
          df_difference > 0 ~
            pchisq(
              chi_square,
              df =
                df_difference,
              lower.tail =
                FALSE
            ),
          
          TRUE ~
            NA_real_
          
        )
        
      )
      
    }
    
  )
  
} else {
  
  likelihood_ratio_tests <- tibble(
    
    reduced_model = character(),
    full_model = character(),
    chi_square = double(),
    df_difference = integer(),
    p_value = double()
    
  )
  
}


# ============================================================
# SELECT A PARSIMONIOUS MODEL
# ============================================================

acceptable_models <- model_comparison |>
  
  filter(
    
    delta_AIC <=
      aic_parsimony_threshold,
    
    convergence ==
      "ok"
    
  )


# Prefer a non-singular model where possible.

non_singular_acceptable_models <-
  acceptable_models |>
  filter(
    !singular
  )


if (
  nrow(
    non_singular_acceptable_models
  ) > 0
) {
  
  acceptable_models <-
    non_singular_acceptable_models
  
}


if (nrow(acceptable_models) == 0) {
  
  acceptable_models <-
    model_comparison |>
    arrange(
      AIC,
      model_order
    ) |>
    slice_head(
      n = 1
    )
  
}


selected_model_name <- acceptable_models |>
  
  arrange(
    model_order
  ) |>
  
  slice_head(
    n = 1
  ) |>
  
  pull(Model)


selected_model_ml <-
  models_ml[[selected_model_name]]


selected_formula <-
  formula(
    selected_model_ml
  )


cat(
  "\nSelected model:",
  selected_model_name,
  "\n"
)

cat(
  "Selected formula:",
  paste(
    deparse(
      selected_formula
    ),
    collapse = ""
  ),
  "\n"
)


# ============================================================
# REFIT SELECTED MODEL USING REML
# ============================================================

final_model <- lmerTest::lmer(
  
  formula =
    selected_formula,
  
  data =
    comparison_data,
  
  weights =
    model_weight,
  
  REML =
    TRUE,
  
  control =
    model_control,
  
  na.action =
    na.fail
  
)


print(
  summary(
    final_model
  )
)


# ============================================================
# SAVE MODEL OBJECTS
# ============================================================

saveRDS(
  models_ml,
  candidate_models_file
)


saveRDS(
  selected_model_ml,
  final_model_ml_file
)


saveRDS(
  final_model,
  final_model_file
)


# ============================================================
# FINAL MODEL COEFFICIENTS
# ============================================================

coefficient_matrix <-
  coef(
    summary(
      final_model
    )
  )


coefficient_table <- tibble(
  
  term =
    rownames(
      coefficient_matrix
    ),
  
  estimate =
    coefficient_matrix[
      ,
      "Estimate"
    ],
  
  standard_error =
    coefficient_matrix[
      ,
      "Std. Error"
    ],
  
  degrees_of_freedom = if (
    "df" %in%
    colnames(
      coefficient_matrix
    )
  ) {
    coefficient_matrix[
      ,
      "df"
    ]
  } else {
    NA_real_
  },
  
  t_value = if (
    "t value" %in%
    colnames(
      coefficient_matrix
    )
  ) {
    coefficient_matrix[
      ,
      "t value"
    ]
  } else {
    NA_real_
  },
  
  p_value = if (
    "Pr(>|t|)" %in%
    colnames(
      coefficient_matrix
    )
  ) {
    coefficient_matrix[
      ,
      "Pr(>|t|)"
    ]
  } else {
    NA_real_
  }
  
)


wald_confidence_intervals <- suppressMessages(
  
  confint(
    
    final_model,
    
    parm =
      "beta_",
    
    method =
      "Wald"
    
  )
  
)


confidence_interval_table <- tibble(
  
  term =
    rownames(
      wald_confidence_intervals
    ),
  
  confidence_low =
    wald_confidence_intervals[
      ,
      1
    ],
  
  confidence_high =
    wald_confidence_intervals[
      ,
      2
    ]
  
)


coefficient_table <- coefficient_table |>
  
  left_join(
    confidence_interval_table,
    by = "term"
  )


# ============================================================
# RANDOM EFFECTS AND VARIANCE COMPONENTS
# ============================================================

random_effect_variance <- as.data.frame(
  VarCorr(
    final_model
  )
) |>
  as_tibble()


event_random_effect_matrix <- ranef(
  final_model
)$publicid


event_random_effects <- tibble(
  
  publicid =
    rownames(
      event_random_effect_matrix
    ),
  
  event_random_intercept =
    event_random_effect_matrix[
      ,
      "(Intercept)"
    ]
  
) |>
  
  arrange(
    event_random_intercept
  )


# ============================================================
# MODEL PERFORMANCE
# ============================================================

model_performance <- tryCatch(
  
  performance::model_performance(
    final_model
  ) |>
    as.data.frame() |>
    as_tibble(),
  
  error = function(e) {
    
    tibble(
      performance_error =
        conditionMessage(e)
    )
    
  }
  
)


collinearity_check <- tryCatch(
  
  performance::check_collinearity(
    final_model
  ) |>
    as.data.frame() |>
    as_tibble(),
  
  error = function(e) {
    
    tibble(
      collinearity_error =
        conditionMessage(e)
    )
    
  }
  
)


singularity_result <- tibble(
  
  selected_model =
    selected_model_name,
  
  singular =
    lme4::isSingular(
      final_model,
      tol = 1e-5
    ),
  
  convergence =
    extract_convergence_message(
      final_model
    ),
  
  residual_standard_deviation =
    sigma(
      final_model
    )
  
)


# ============================================================
# PREDICTIONS AND RESIDUALS
# ============================================================

model_sample <- comparison_data |>
  
  mutate(
    
    selected_model =
      selected_model_name,
    
    predicted_population_mmi =
      as.numeric(
        predict(
          final_model,
          re.form = NA
        )
      ),
    
    predicted_event_mmi =
      as.numeric(
        predict(
          final_model,
          re.form = NULL
        )
      ),
    
    residual_mmi =
      reported_mmi -
      predicted_event_mmi,
    
    standardised_residual =
      residual_mmi /
      sigma(
        final_model
      )
    
  )


prediction_fields <- model_sample |>
  
  select(
    
    observation_id,
    
    selected_model,
    
    predicted_population_mmi,
    predicted_event_mmi,
    
    residual_mmi,
    standardised_residual,
    
    model_weight,
    capped_report_count
    
  )


iar_with_model_results <- iar_source |>
  
  left_join(
    prediction_fields,
    by = "observation_id"
  ) |>
  
  mutate(
    
    model_included =
      !is.na(
        predicted_event_mmi
      )
    
  )


# ============================================================
# MODEL EXCLUSION SUMMARY
# ============================================================

comparison_observation_ids <-
  comparison_data$observation_id


model_exclusion_details <- iar_prepared |>
  
  mutate(
    
    model_included =
      observation_id %in%
      comparison_observation_ids,
    
    final_exclusion_reason = case_when(
      
      model_included ~
        "included",
      
      !is.na(exclusion_reason) ~
        exclusion_reason,
      
      use_tectonic_region &
        is.na(
          tectonic_region_model
        ) ~
        "missing_or_sparse_tectonic_region",
      
      use_fault_style &
        is.na(
          fault_style_model
        ) ~
        "missing_or_sparse_fault_style",
      
      TRUE ~
        "excluded_from_common_comparison_sample"
      
    )
    
  )


model_exclusion_summary <- model_exclusion_details |>
  
  count(
    final_exclusion_reason,
    name = "observation_count"
  ) |>
  
  mutate(
    
    percentage =
      100 *
      observation_count /
      sum(
        observation_count
      )
    
  ) |>
  
  arrange(
    desc(
      observation_count
    )
  )


# ============================================================
# OPTIONAL TERM AVAILABILITY SUMMARY
# ============================================================

optional_term_summary <- tibble(
  
  predictor = c(
    "tectonic_region",
    "fault_style"
  ),
  
  requested = c(
    request_tectonic_region_model,
    request_fault_style_model
  ),
  
  initially_usable = c(
    tectonic_assessment$usable,
    fault_assessment$usable
  ),
  
  retained_fraction = c(
    tectonic_assessment$retained_fraction,
    fault_assessment$retained_fraction
  ),
  
  retained_level_count = c(
    length(
      tectonic_assessment$retained_levels
    ),
    length(
      fault_assessment$retained_levels
    )
  ),
  
  used_in_models = c(
    use_tectonic_region,
    use_fault_style
  )
  
)


# ============================================================
# SAVE TABLE OUTPUTS
# ============================================================

arrow::write_parquet(
  model_sample,
  model_sample_file
)


arrow::write_parquet(
  iar_with_model_results,
  model_results_file
)


readr::write_csv(
  model_comparison,
  file.path(
    result_directory,
    "model_comparison.csv"
  )
)


readr::write_csv(
  likelihood_ratio_tests,
  file.path(
    result_directory,
    "model_likelihood_ratio_tests.csv"
  )
)


readr::write_csv(
  coefficient_table,
  file.path(
    result_directory,
    "final_model_coefficients.csv"
  )
)


readr::write_csv(
  random_effect_variance,
  file.path(
    result_directory,
    "final_model_variance_components.csv"
  )
)


readr::write_csv(
  event_random_effects,
  file.path(
    result_directory,
    "event_random_effects.csv"
  )
)


readr::write_csv(
  model_performance,
  file.path(
    result_directory,
    "final_model_performance.csv"
  )
)


readr::write_csv(
  collinearity_check,
  file.path(
    result_directory,
    "final_model_collinearity.csv"
  )
)


readr::write_csv(
  singularity_result,
  file.path(
    result_directory,
    "final_model_fit_status.csv"
  )
)


readr::write_csv(
  optional_term_summary,
  file.path(
    result_directory,
    "optional_predictor_availability.csv"
  )
)


readr::write_csv(
  tectonic_assessment$summary,
  file.path(
    result_directory,
    "tectonic_region_level_summary.csv"
  )
)


readr::write_csv(
  fault_assessment$summary,
  file.path(
    result_directory,
    "fault_style_level_summary.csv"
  )
)


readr::write_csv(
  model_exclusion_summary,
  file.path(
    diagnostic_directory,
    "model_exclusion_summary.csv"
  )
)


capture.output(
  
  summary(
    final_model
  ),
  
  file =
    file.path(
      result_directory,
      "final_model_summary.txt"
    )
  
)


# ============================================================
# DIAGNOSTIC PLOTS
# ============================================================

residual_distance_plot <- ggplot(
  
  model_sample,
  
  aes(
    x = model_distance_km,
    y = residual_mmi
  )
  
) +
  
  geom_point(
    alpha = 0.15
  ) +
  
  geom_hline(
    yintercept = 0
  ) +
  
  scale_x_log10() +
  
  theme_bw() +
  
  labs(
    
    title =
      "MMI residuals against model distance",
    
    subtitle =
      paste(
        "Distance variable:",
        distance_variable
      ),
    
    x =
      paste0(
        distance_variable,
        " (km)"
      ),
    
    y =
      "Conditional MMI residual"
    
  )


residual_vs30_plot <- ggplot(
  
  model_sample,
  
  aes(
    x = Vs30_m_s,
    y = residual_mmi
  )
  
) +
  
  geom_point(
    alpha = 0.15
  ) +
  
  geom_hline(
    yintercept = 0
  ) +
  
  scale_x_log10() +
  
  theme_bw() +
  
  labs(
    
    title =
      "MMI residuals against Vs30",
    
    x =
      "Vs30 (m/s)",
    
    y =
      "Conditional MMI residual"
    
  )


residual_fitted_plot <- ggplot(
  
  model_sample,
  
  aes(
    x = predicted_event_mmi,
    y = residual_mmi
  )
  
) +
  
  geom_point(
    alpha = 0.15
  ) +
  
  geom_hline(
    yintercept = 0
  ) +
  
  theme_bw() +
  
  labs(
    
    title =
      "Residuals against fitted MMI",
    
    x =
      "Conditional fitted MMI",
    
    y =
      "Conditional MMI residual"
    
  )


observed_predicted_plot <- ggplot(
  
  model_sample,
  
  aes(
    x = predicted_event_mmi,
    y = reported_mmi
  )
  
) +
  
  geom_point(
    alpha = 0.15
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0
  ) +
  
  coord_equal() +
  
  theme_bw() +
  
  labs(
    
    title =
      "Observed versus predicted MMI",
    
    x =
      "Conditional predicted MMI",
    
    y =
      "Observed Felt RAPID MMI"
    
  )


residual_qq_plot <- ggplot(
  
  model_sample,
  
  aes(
    sample =
      standardised_residual
  )
  
) +
  
  stat_qq(
    alpha = 0.25
  ) +
  
  stat_qq_line() +
  
  theme_bw() +
  
  labs(
    
    title =
      "Normal Q-Q plot of standardised residuals",
    
    x =
      "Theoretical quantile",
    
    y =
      "Observed standardised residual"
    
  )


event_effect_plot <- event_random_effects |>
  
  mutate(
    
    publicid =
      forcats::fct_reorder(
        publicid,
        event_random_intercept
      )
    
  ) |>
  
  ggplot(
    
    aes(
      x = publicid,
      y = event_random_intercept
    )
    
  ) +
  
  geom_point() +
  
  geom_hline(
    yintercept = 0
  ) +
  
  coord_flip() +
  
  theme_bw() +
  
  theme(
    axis.text.y =
      element_blank(),
    axis.ticks.y =
      element_blank()
  ) +
  
  labs(
    
    title =
      "Earthquake random intercepts",
    
    x =
      "Earthquake",
    
    y =
      "Event MMI adjustment"
    
  )


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "residuals_vs_distance.png"
    ),
  
  plot =
    residual_distance_plot,
  
  width = 8,
  height = 6,
  dpi = 300
  
)


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "residuals_vs_vs30.png"
    ),
  
  plot =
    residual_vs30_plot,
  
  width = 8,
  height = 6,
  dpi = 300
  
)


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "residuals_vs_fitted.png"
    ),
  
  plot =
    residual_fitted_plot,
  
  width = 8,
  height = 6,
  dpi = 300
  
)


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "observed_vs_predicted.png"
    ),
  
  plot =
    observed_predicted_plot,
  
  width = 7,
  height = 7,
  dpi = 300
  
)


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "residual_qq_plot.png"
    ),
  
  plot =
    residual_qq_plot,
  
  width = 7,
  height = 7,
  dpi = 300
  
)


ggsave(
  
  filename =
    file.path(
      plot_directory,
      "event_random_effects.png"
    ),
  
  plot =
    event_effect_plot,
  
  width = 8,
  height = 10,
  dpi = 300,
  limitsize = FALSE
  
)


# ============================================================
# PERFORMANCE CHECK OBJECT
# ============================================================

model_check <- tryCatch(
  
  performance::check_model(
    
    final_model,
    
    show_dots = FALSE
    
  ),
  
  error = function(e) {
    e
  }
  
)


saveRDS(
  
  model_check,
  
  file.path(
    result_directory,
    "final_model_check.rds"
  )
  
)


# ============================================================
# VERIFY IMPORTANT OUTPUT FILES
# ============================================================

important_output_files <- c(
  
  final_model_file,
  final_model_ml_file,
  candidate_models_file,
  
  model_sample_file,
  model_results_file,
  
  file.path(
    result_directory,
    "model_comparison.csv"
  ),
  
  file.path(
    result_directory,
    "final_model_coefficients.csv"
  ),
  
  file.path(
    result_directory,
    "final_model_performance.csv"
  ),
  
  file.path(
    result_directory,
    "final_model_summary.txt"
  ),
  
  file.path(
    plot_directory,
    "residuals_vs_distance.png"
  ),
  
  file.path(
    plot_directory,
    "observed_vs_predicted.png"
  )
  
)


invalid_output_files <- important_output_files[
  
  !file.exists(
    important_output_files
  ) |
    
    is.na(
      file.info(
        important_output_files
      )$size
    ) |
    
    file.info(
      important_output_files
    )$size == 0
  
]


if (length(invalid_output_files) > 0) {
  
  stop(
    paste0(
      "One or more Script 8 outputs were not written correctly:\n",
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
  "\nNZ Felt RAPID attenuation modelling complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Input observations:",
  format(
    nrow(iar_source),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Model observations:",
  format(
    nrow(model_sample),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Model earthquakes:",
  format(
    n_distinct(
      model_sample$publicid
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Selected model:",
  selected_model_name,
  "\n"
)

cat(
  "Distance variable:",
  distance_variable,
  "\n"
)

cat(
  "Model is singular:",
  lme4::isSingular(
    final_model,
    tol = 1e-5
  ),
  "\n"
)

cat(
  "Convergence status:",
  extract_convergence_message(
    final_model
  ),
  "\n"
)

cat(
  "Residual standard deviation:",
  round(
    sigma(
      final_model
    ),
    4
  ),
  "\n"
)

cat(
  "\nFinal REML model:\n",
  normalizePath(
    final_model_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nModel sample:\n",
  normalizePath(
    model_sample_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nDataset with model predictions and residuals:\n",
  normalizePath(
    model_results_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nModel result tables:\n",
  normalizePath(
    result_directory,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)