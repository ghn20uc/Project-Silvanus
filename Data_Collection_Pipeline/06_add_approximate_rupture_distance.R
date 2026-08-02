# ============================================================
# 06_add_approximate_rupture_distance.R
#
# Estimate finite-rupture distance for the NZ Felt RAPID
# attenuation dataset.
#
# This script calculates candidate rupture distances for both
# CMT nodal planes. These are approximate distances because:
#
#   1. rupture dimensions are predicted from Mw;
#   2. the rectangular rupture is centred on the event/CMT
#      centre;
#   3. the physical rupture plane is not known from the
#      moment tensor alone;
#   4. mapped finite-fault models are not yet incorporated.
#
# Required input:
#   data_processed/NZ_FR_IAR_dataset.parquet
#
# Outputs:
#   data_processed/rupture_geometry.parquet
#   data_processed/NZ_FR_IAR_with_Rrup.parquet
#   data_processed/diagnostics/iar_rrup_*.csv
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(tidyverse)
library(arrow)
library(geosphere)


# ============================================================
# RUPTURE-DISTANCE CONTROLS
# ============================================================

# Script 6 makes no API requests, so date and magnitude download
# controls are not needed here.
#
# The following settings control the rupture approximation.


# Use fault-style-specific Wells-Coppersmith coefficients where
# a reliable broad fault style is available.
#
# FALSE:
#   use the all-slip coefficients for every eligible event.
#
# TRUE:
#   use strike-slip, reverse, or normal coefficients where
#   available, and otherwise use the all-slip coefficients.
#
# FALSE is the safer starting point because many events do not
# have an unambiguous physical rupture plane.

use_fault_style_specific_scaling <- FALSE


# Prevent use of the scaling relationships outside their
# approximate calibration magnitude ranges.
#
# FALSE is recommended.

allow_scaling_extrapolation <- FALSE


# Rupture centre:
#
# "cmt_then_hypocentre"
#   use CMT longitude, latitude and centroid depth where all
#   three are available; otherwise use catalogue hypocentre.
#
# "hypocentre"
#   always use catalogue event coordinates and depth.

rupture_centre_preference <- "cmt_then_hypocentre"


# Candidate used for sensitivity-analysis output:
#
# "minimum" = smaller of the two nodal-plane candidates
# "mean"    = arithmetic mean of available candidates
# "np1"     = nodal plane 1
# "np2"     = nodal plane 2
#
# This does not discard either individual candidate.

rrup_candidate_for_sensitivity <- "minimum"


# Do not replace the main Rhypo-based model distance by default.
# Set TRUE only for a deliberate finite-rupture sensitivity
# model.

use_approx_rrup_as_primary_distance <- FALSE


# ------------------------------------------------------------
# Validate controls
# ------------------------------------------------------------

valid_centre_preferences <- c(
  "cmt_then_hypocentre",
  "hypocentre"
)


if (
  !rupture_centre_preference %in%
  valid_centre_preferences
) {
  
  stop(
    paste0(
      "rupture_centre_preference must be one of:\n",
      paste(
        valid_centre_preferences,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


valid_rrup_candidates <- c(
  "minimum",
  "mean",
  "np1",
  "np2"
)


if (
  !rrup_candidate_for_sensitivity %in%
  valid_rrup_candidates
) {
  
  stop(
    paste0(
      "rrup_candidate_for_sensitivity must be one of:\n",
      paste(
        valid_rrup_candidates,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# ============================================================
# PATHS
# ============================================================

input_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_dataset.parquet"
)

rupture_geometry_file <- file.path(
  "data_processed",
  "rupture_geometry.parquet"
)

output_file <- file.path(
  "data_processed",
  "NZ_FR_IAR_with_Rrup.parquet"
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


# ============================================================
# LOAD SCRIPT 5 DATASET
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Script 5 output was not found:\n",
      normalizePath(
        input_file,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun Script 5 first."
    ),
    call. = FALSE
  )
  
}


iar <- arrow::read_parquet(
  input_file
) |>
  as_tibble()


if (nrow(iar) == 0) {
  
  stop(
    "The Script 5 modelling dataset contains zero rows.",
    call. = FALSE
  )
  
}


cat(
  "Felt RAPID modelling rows loaded:",
  format(
    nrow(iar),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "publicid",
  "report_longitude",
  "report_latitude",
  "event_longitude",
  "event_latitude",
  "event_depth_km",
  "Repi_km",
  "Rhypo_km"
)


missing_columns <- setdiff(
  required_columns,
  names(iar)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "The Script 5 dataset is missing required columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


# ============================================================
# ADD OPTIONAL COLUMNS WHERE ABSENT
# ============================================================

optional_numeric_columns <- c(
  "Mw_model",
  "mt_longitude",
  "mt_latitude",
  "mt_centroid_depth_km",
  "mt_strike1",
  "mt_dip1",
  "mt_rake1",
  "mt_strike2",
  "mt_dip2",
  "mt_rake2",
  "R_model_km"
)


for (column_name in optional_numeric_columns) {
  
  if (!column_name %in% names(iar)) {
    
    iar[[column_name]] <- NA_real_
    
  }
  
}


optional_character_columns <- c(
  "mt_fault_style_consensus",
  "fault_style",
  "R_model_type"
)


for (column_name in optional_character_columns) {
  
  if (!column_name %in% names(iar)) {
    
    iar[[column_name]] <- NA_character_
    
  }
  
}


# ============================================================
# STANDARDISE TYPES
# ============================================================

numeric_columns_to_parse <- intersect(
  
  c(
    "report_longitude",
    "report_latitude",
    "event_longitude",
    "event_latitude",
    "event_depth_km",
    "Repi_km",
    "Rhypo_km",
    "Mw_model",
    "mt_longitude",
    "mt_latitude",
    "mt_centroid_depth_km",
    "mt_strike1",
    "mt_dip1",
    "mt_rake1",
    "mt_strike2",
    "mt_dip2",
    "mt_rake2",
    "R_model_km"
  ),
  
  names(iar)
  
)


iar <- iar |>
  
  mutate(
    
    publicid = publicid |>
      as.character() |>
      stringr::str_trim() |>
      na_if(""),
    
    across(
      all_of(numeric_columns_to_parse),
      ~ suppressWarnings(
        as.numeric(.x)
      )
    )
    
  )


# ============================================================
# BUILD ONE-ROW-PER-EARTHQUAKE RUPTURE TABLE
# ============================================================

rupture_events <- iar |>
  
  arrange(
    publicid
  ) |>
  
  distinct(
    publicid,
    .keep_all = TRUE
  ) |>
  
  transmute(
    
    publicid,
    
    rupture_mw = Mw_model,
    
    event_longitude,
    event_latitude,
    event_depth_km,
    
    mt_longitude,
    mt_latitude,
    mt_centroid_depth_km,
    
    strike_np1_deg = mt_strike1,
    dip_np1_deg = mt_dip1,
    rake_np1_deg = mt_rake1,
    
    strike_np2_deg = mt_strike2,
    dip_np2_deg = mt_dip2,
    rake_np2_deg = mt_rake2,
    
    rupture_fault_style_raw = coalesce(
      mt_fault_style_consensus,
      fault_style
    )
    
  ) |>
  
  mutate(
    
    rupture_fault_style_raw =
      rupture_fault_style_raw |>
      as.character() |>
      stringr::str_to_lower() |>
      stringr::str_trim(),
    
    rupture_fault_style = case_when(
      
      stringr::str_detect(
        rupture_fault_style_raw,
        "strike"
      ) ~
        "strike-slip",
      
      stringr::str_detect(
        rupture_fault_style_raw,
        "reverse|thrust"
      ) ~
        "reverse",
      
      stringr::str_detect(
        rupture_fault_style_raw,
        "normal"
      ) ~
        "normal",
      
      TRUE ~
        "unknown"
      
    ),
    
    valid_mt_centre =
      !is.na(mt_longitude) &
      !is.na(mt_latitude) &
      !is.na(mt_centroid_depth_km) &
      mt_longitude >= -180 &
      mt_longitude <= 180 &
      mt_latitude >= -90 &
      mt_latitude <= 90 &
      mt_centroid_depth_km >= 0 &
      mt_centroid_depth_km <= 700,
    
    valid_hypocentre =
      !is.na(event_longitude) &
      !is.na(event_latitude) &
      !is.na(event_depth_km) &
      event_longitude >= -180 &
      event_longitude <= 180 &
      event_latitude >= -90 &
      event_latitude <= 90 &
      event_depth_km >= 0 &
      event_depth_km <= 700,
    
    .use_mt_centre =
      rupture_centre_preference ==
      "cmt_then_hypocentre" &
      valid_mt_centre,
    
    rupture_centre_longitude = case_when(
      
      .use_mt_centre ~
        mt_longitude,
      
      valid_hypocentre ~
        event_longitude,
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_centre_latitude = case_when(
      
      .use_mt_centre ~
        mt_latitude,
      
      valid_hypocentre ~
        event_latitude,
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_centre_depth_km = case_when(
      
      .use_mt_centre ~
        mt_centroid_depth_km,
      
      valid_hypocentre ~
        event_depth_km,
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_centre_source = case_when(
      
      .use_mt_centre ~
        "GeoNet_CMT_centroid",
      
      valid_hypocentre ~
        "GeoNet_catalogue_hypocentre",
      
      TRUE ~
        NA_character_
      
    )
    
  ) |>
  
  select(
    -.use_mt_centre
  )


# ============================================================
# ASSIGN WELLS-COPPERSMITH SCALING COEFFICIENTS
#
# RLD = subsurface rupture length
# RW  = downdip rupture width
#
# log10(RLD) = a_length + b_length * Mw
# log10(RW)  = a_width  + b_width  * Mw
# ============================================================

if (use_fault_style_specific_scaling) {
  
  rupture_events <- rupture_events |>
    
    mutate(
      
      scaling_relation = case_when(
        
        rupture_fault_style ==
          "strike-slip" ~
          "strike-slip",
        
        rupture_fault_style ==
          "reverse" ~
          "reverse",
        
        rupture_fault_style ==
          "normal" ~
          "normal",
        
        TRUE ~
          "all-slip"
        
      ),
      
      length_coefficient_a = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          -2.57,
        
        scaling_relation ==
          "reverse" ~
          -2.42,
        
        scaling_relation ==
          "normal" ~
          -1.88,
        
        TRUE ~
          -2.44
        
      ),
      
      length_coefficient_b = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          0.62,
        
        scaling_relation ==
          "reverse" ~
          0.58,
        
        scaling_relation ==
          "normal" ~
          0.50,
        
        TRUE ~
          0.59
        
      ),
      
      width_coefficient_a = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          -0.76,
        
        scaling_relation ==
          "reverse" ~
          -1.61,
        
        scaling_relation ==
          "normal" ~
          -1.14,
        
        TRUE ~
          -1.01
        
      ),
      
      width_coefficient_b = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          0.27,
        
        scaling_relation ==
          "reverse" ~
          0.41,
        
        scaling_relation ==
          "normal" ~
          0.35,
        
        TRUE ~
          0.32
        
      ),
      
      scaling_mw_min = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          4.8,
        
        scaling_relation ==
          "reverse" ~
          4.8,
        
        scaling_relation ==
          "normal" ~
          5.2,
        
        TRUE ~
          4.8
        
      ),
      
      scaling_mw_max = case_when(
        
        scaling_relation ==
          "strike-slip" ~
          8.1,
        
        scaling_relation ==
          "reverse" ~
          7.6,
        
        scaling_relation ==
          "normal" ~
          7.3,
        
        TRUE ~
          8.1
        
      )
      
    )
  
} else {
  
  rupture_events <- rupture_events |>
    
    mutate(
      
      scaling_relation =
        "all-slip",
      
      length_coefficient_a =
        -2.44,
      
      length_coefficient_b =
        0.59,
      
      width_coefficient_a =
        -1.01,
      
      width_coefficient_b =
        0.32,
      
      scaling_mw_min =
        4.8,
      
      scaling_mw_max =
        8.1
      
    )
  
}


# ============================================================
# CALCULATE RUPTURE DIMENSIONS
# ============================================================

rupture_events <- rupture_events |>
  
  mutate(
    
    rupture_scaling_extrapolated =
      !is.na(rupture_mw) &
      (
        rupture_mw < scaling_mw_min |
          rupture_mw > scaling_mw_max
      ),
    
    rupture_scaling_eligible =
      !is.na(rupture_mw) &
      rupture_mw > 0 &
      (
        allow_scaling_extrapolation |
          !rupture_scaling_extrapolated
      ),
    
    rupture_length_km = case_when(
      
      rupture_scaling_eligible ~
        10^(
          length_coefficient_a +
            length_coefficient_b *
            rupture_mw
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    predicted_rupture_width_km = case_when(
      
      rupture_scaling_eligible ~
        10^(
          width_coefficient_a +
            width_coefficient_b *
            rupture_mw
        ),
      
      TRUE ~
        NA_real_
      
    )
    
  )


# ============================================================
# CAP WIDTH SO A CENTRED RECTANGLE DOES NOT EXTEND ABOVE
# THE GROUND SURFACE
# ============================================================

rupture_events <- rupture_events |>
  
  mutate(
    
    valid_np1_inputs =
      rupture_scaling_eligible &
      !is.na(rupture_centre_longitude) &
      !is.na(rupture_centre_latitude) &
      !is.na(rupture_centre_depth_km) &
      !is.na(strike_np1_deg) &
      strike_np1_deg >= 0 &
      strike_np1_deg <= 360 &
      !is.na(dip_np1_deg) &
      dip_np1_deg > 0 &
      dip_np1_deg <= 90 &
      !is.na(rupture_length_km) &
      rupture_length_km > 0 &
      !is.na(predicted_rupture_width_km) &
      predicted_rupture_width_km > 0,
    
    valid_np2_inputs =
      rupture_scaling_eligible &
      !is.na(rupture_centre_longitude) &
      !is.na(rupture_centre_latitude) &
      !is.na(rupture_centre_depth_km) &
      !is.na(strike_np2_deg) &
      strike_np2_deg >= 0 &
      strike_np2_deg <= 360 &
      !is.na(dip_np2_deg) &
      dip_np2_deg > 0 &
      dip_np2_deg <= 90 &
      !is.na(rupture_length_km) &
      rupture_length_km > 0 &
      !is.na(predicted_rupture_width_km) &
      predicted_rupture_width_km > 0,
    
    maximum_centred_width_np1_km = case_when(
      
      valid_np1_inputs ~
        2 *
        rupture_centre_depth_km /
        sin(
          dip_np1_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    maximum_centred_width_np2_km = case_when(
      
      valid_np2_inputs ~
        2 *
        rupture_centre_depth_km /
        sin(
          dip_np2_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    effective_rupture_width_np1_km = case_when(
      
      valid_np1_inputs ~
        pmin(
          predicted_rupture_width_km,
          maximum_centred_width_np1_km
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    effective_rupture_width_np2_km = case_when(
      
      valid_np2_inputs ~
        pmin(
          predicted_rupture_width_km,
          maximum_centred_width_np2_km
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_width_capped_np1 =
      valid_np1_inputs &
      effective_rupture_width_np1_km <
      predicted_rupture_width_km,
    
    rupture_width_capped_np2 =
      valid_np2_inputs &
      effective_rupture_width_np2_km <
      predicted_rupture_width_km,
    
    rupture_geometry_available_np1 =
      valid_np1_inputs &
      !is.na(effective_rupture_width_np1_km) &
      effective_rupture_width_np1_km > 0,
    
    rupture_geometry_available_np2 =
      valid_np2_inputs &
      !is.na(effective_rupture_width_np2_km) &
      effective_rupture_width_np2_km > 0,
    
    rupture_top_depth_np1_km = case_when(
      
      rupture_geometry_available_np1 ~
        rupture_centre_depth_km -
        0.5 *
        effective_rupture_width_np1_km *
        sin(
          dip_np1_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_bottom_depth_np1_km = case_when(
      
      rupture_geometry_available_np1 ~
        rupture_centre_depth_km +
        0.5 *
        effective_rupture_width_np1_km *
        sin(
          dip_np1_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_top_depth_np2_km = case_when(
      
      rupture_geometry_available_np2 ~
        rupture_centre_depth_km -
        0.5 *
        effective_rupture_width_np2_km *
        sin(
          dip_np2_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    ),
    
    rupture_bottom_depth_np2_km = case_when(
      
      rupture_geometry_available_np2 ~
        rupture_centre_depth_km +
        0.5 *
        effective_rupture_width_np2_km *
        sin(
          dip_np2_deg *
            pi /
            180
        ),
      
      TRUE ~
        NA_real_
      
    )
    
  )


# ============================================================
# SAVE EVENT-LEVEL RUPTURE GEOMETRY
# ============================================================

arrow::write_parquet(
  rupture_events,
  rupture_geometry_file
)


# ============================================================
# JOIN RUPTURE GEOMETRY TO FELT RAPID OBSERVATIONS
# ============================================================

rupture_columns_for_join <- rupture_events |>
  
  select(
    publicid,
    rupture_mw,
    rupture_fault_style,
    scaling_relation,
    scaling_mw_min,
    scaling_mw_max,
    rupture_scaling_extrapolated,
    rupture_scaling_eligible,
    rupture_centre_longitude,
    rupture_centre_latitude,
    rupture_centre_depth_km,
    rupture_centre_source,
    rupture_length_km,
    predicted_rupture_width_km,
    strike_np1_deg,
    dip_np1_deg,
    rake_np1_deg,
    effective_rupture_width_np1_km,
    rupture_width_capped_np1,
    rupture_top_depth_np1_km,
    rupture_bottom_depth_np1_km,
    rupture_geometry_available_np1,
    strike_np2_deg,
    dip_np2_deg,
    rake_np2_deg,
    effective_rupture_width_np2_km,
    rupture_width_capped_np2,
    rupture_top_depth_np2_km,
    rupture_bottom_depth_np2_km,
    rupture_geometry_available_np2
  )


iar_rrup <- iar |>
  
  left_join(
    rupture_columns_for_join,
    by = "publicid"
  )


if (nrow(iar_rrup) != nrow(iar)) {
  
  stop(
    paste0(
      "The rupture-geometry join changed the number of rows.\n",
      "Before join: ",
      nrow(iar),
      "\nAfter join: ",
      nrow(iar_rrup)
    ),
    call. = FALSE
  )
  
}


# ============================================================
# CALCULATE CENTRE-TO-SITE HORIZONTAL COMPONENTS
#
# Bearing is measured clockwise from north.
# x is positive east and y is positive north.
# ============================================================

iar_rrup$rupture_centre_Repi_km <-
  NA_real_

iar_rrup$rupture_centre_bearing_deg <-
  NA_real_

iar_rrup$site_east_of_centre_km <-
  NA_real_

iar_rrup$site_north_of_centre_km <-
  NA_real_


valid_horizontal_rows <- which(
  
  !is.na(iar_rrup$rupture_centre_longitude) &
    !is.na(iar_rrup$rupture_centre_latitude) &
    !is.na(iar_rrup$report_longitude) &
    !is.na(iar_rrup$report_latitude) &
    
    iar_rrup$rupture_centre_longitude >= -180 &
    iar_rrup$rupture_centre_longitude <= 180 &
    
    iar_rrup$rupture_centre_latitude >= -90 &
    iar_rrup$rupture_centre_latitude <= 90 &
    
    iar_rrup$report_longitude >= -180 &
    iar_rrup$report_longitude <= 180 &
    
    iar_rrup$report_latitude >= -90 &
    iar_rrup$report_latitude <= 90
  
)


if (length(valid_horizontal_rows) > 0) {
  
  centre_coordinates <- cbind(
    
    iar_rrup$rupture_centre_longitude[
      valid_horizontal_rows
    ],
    
    iar_rrup$rupture_centre_latitude[
      valid_horizontal_rows
    ]
    
  )
  
  
  report_coordinates <- cbind(
    
    iar_rrup$report_longitude[
      valid_horizontal_rows
    ],
    
    iar_rrup$report_latitude[
      valid_horizontal_rows
    ]
    
  )
  
  
  centre_distance_km <-
    geosphere::distGeo(
      centre_coordinates,
      report_coordinates
    ) /
    1000
  
  
  centre_bearing_deg <-
    geosphere::bearing(
      centre_coordinates,
      report_coordinates
    )
  
  
  centre_bearing_deg[
    is.na(centre_bearing_deg) &
      centre_distance_km == 0
  ] <- 0
  
  
  centre_bearing_rad <-
    centre_bearing_deg *
    pi /
    180
  
  
  iar_rrup$rupture_centre_Repi_km[
    valid_horizontal_rows
  ] <- centre_distance_km
  
  
  iar_rrup$rupture_centre_bearing_deg[
    valid_horizontal_rows
  ] <- centre_bearing_deg
  
  
  iar_rrup$site_east_of_centre_km[
    valid_horizontal_rows
  ] <-
    centre_distance_km *
    sin(
      centre_bearing_rad
    )
  
  
  iar_rrup$site_north_of_centre_km[
    valid_horizontal_rows
  ] <-
    centre_distance_km *
    cos(
      centre_bearing_rad
    )
  
}


# ============================================================
# DISTANCE FROM A SITE TO A FINITE RECTANGULAR PLANE
#
# Coordinate system:
#   x = east
#   y = north
#   z = positive downward
#
# The rupture rectangle is centred at the selected rupture
# centre.
# ============================================================

distance_to_rupture_rectangle <- function(
    site_east_km,
    site_north_km,
    centre_depth_km,
    rupture_length_km,
    rupture_width_km,
    strike_deg,
    dip_deg
) {
  
  output <- rep(
    NA_real_,
    length(site_east_km)
  )
  
  
  valid <- is.finite(site_east_km) &
    is.finite(site_north_km) &
    is.finite(centre_depth_km) &
    centre_depth_km >= 0 &
    is.finite(rupture_length_km) &
    rupture_length_km > 0 &
    is.finite(rupture_width_km) &
    rupture_width_km > 0 &
    is.finite(strike_deg) &
    strike_deg >= 0 &
    strike_deg <= 360 &
    is.finite(dip_deg) &
    dip_deg > 0 &
    dip_deg <= 90
  
  
  if (!any(valid)) {
    
    return(output)
    
  }
  
  
  strike_rad <-
    strike_deg[valid] *
    pi /
    180
  
  
  dip_rad <-
    dip_deg[valid] *
    pi /
    180
  
  
  # Unit vector along strike
  strike_east <-
    sin(strike_rad)
  
  strike_north <-
    cos(strike_rad)
  
  strike_down <-
    rep(
      0,
      sum(valid)
    )
  
  
  # Unit vector down dip, to the right of strike
  dip_east <-
    cos(strike_rad) *
    cos(dip_rad)
  
  dip_north <-
    -sin(strike_rad) *
    cos(dip_rad)
  
  dip_down <-
    sin(dip_rad)
  
  
  # Vector from rupture centre to the surface site
  site_vector_east <-
    site_east_km[valid]
  
  site_vector_north <-
    site_north_km[valid]
  
  site_vector_down <-
    -centre_depth_km[valid]
  
  
  # Coordinates of the site projection in the fault-plane
  # basis
  along_strike_projection <-
    site_vector_east *
    strike_east +
    site_vector_north *
    strike_north +
    site_vector_down *
    strike_down
  
  
  down_dip_projection <-
    site_vector_east *
    dip_east +
    site_vector_north *
    dip_north +
    site_vector_down *
    dip_down
  
  
  half_length <-
    rupture_length_km[valid] /
    2
  
  
  half_width <-
    rupture_width_km[valid] /
    2
  
  
  # Restrict the closest point to the finite rectangle
  closest_along_strike <-
    pmin(
      pmax(
        along_strike_projection,
        -half_length
      ),
      half_length
    )
  
  
  closest_down_dip <-
    pmin(
      pmax(
        down_dip_projection,
        -half_width
      ),
      half_width
    )
  
  
  residual_east <-
    site_vector_east -
    closest_along_strike *
    strike_east -
    closest_down_dip *
    dip_east
  
  
  residual_north <-
    site_vector_north -
    closest_along_strike *
    strike_north -
    closest_down_dip *
    dip_north
  
  
  residual_down <-
    site_vector_down -
    closest_along_strike *
    strike_down -
    closest_down_dip *
    dip_down
  
  
  output[valid] <- sqrt(
    
    residual_east^2 +
      residual_north^2 +
      residual_down^2
    
  )
  
  
  output
  
}


# ============================================================
# CALCULATE CANDIDATE RRUP FOR BOTH NODAL PLANES
# ============================================================

iar_rrup <- iar_rrup |>
  
  mutate(
    
    Rrup_approx_np1_km =
      distance_to_rupture_rectangle(
        
        site_east_km =
          site_east_of_centre_km,
        
        site_north_km =
          site_north_of_centre_km,
        
        centre_depth_km =
          rupture_centre_depth_km,
        
        rupture_length_km =
          rupture_length_km,
        
        rupture_width_km =
          effective_rupture_width_np1_km,
        
        strike_deg =
          strike_np1_deg,
        
        dip_deg =
          dip_np1_deg
        
      ),
    
    Rrup_approx_np2_km =
      distance_to_rupture_rectangle(
        
        site_east_km =
          site_east_of_centre_km,
        
        site_north_km =
          site_north_of_centre_km,
        
        centre_depth_km =
          rupture_centre_depth_km,
        
        rupture_length_km =
          rupture_length_km,
        
        rupture_width_km =
          effective_rupture_width_np2_km,
        
        strike_deg =
          strike_np2_deg,
        
        dip_deg =
          dip_np2_deg
        
      )
    
  )


# ============================================================
# COMBINE CANDIDATE RESULTS WITHOUT LOSING PLANE AMBIGUITY
# ============================================================

pair_min_na <- function(x, y) {
  
  case_when(
    
    is.na(x) &
      is.na(y) ~
      NA_real_,
    
    is.na(x) ~
      y,
    
    is.na(y) ~
      x,
    
    TRUE ~
      pmin(
        x,
        y
      )
    
  )
  
}


pair_max_na <- function(x, y) {
  
  case_when(
    
    is.na(x) &
      is.na(y) ~
      NA_real_,
    
    is.na(x) ~
      y,
    
    is.na(y) ~
      x,
    
    TRUE ~
      pmax(
        x,
        y
      )
    
  )
  
}


pair_mean_na <- function(x, y) {
  
  case_when(
    
    is.na(x) &
      is.na(y) ~
      NA_real_,
    
    is.na(x) ~
      y,
    
    is.na(y) ~
      x,
    
    TRUE ~
      (
        x +
          y
      ) /
      2
    
  )
  
}


iar_rrup <- iar_rrup |>
  
  mutate(
    
    Rrup_approx_min_km =
      pair_min_na(
        Rrup_approx_np1_km,
        Rrup_approx_np2_km
      ),
    
    Rrup_approx_max_km =
      pair_max_na(
        Rrup_approx_np1_km,
        Rrup_approx_np2_km
      ),
    
    Rrup_approx_mean_km =
      pair_mean_na(
        Rrup_approx_np1_km,
        Rrup_approx_np2_km
      ),
    
    Rrup_nodal_plane_range_km =
      Rrup_approx_max_km -
      Rrup_approx_min_km,
    
    Rrup_approx_plane_basis = case_when(
      
      !is.na(Rrup_approx_np1_km) &
        !is.na(Rrup_approx_np2_km) ~
        "both_candidate_nodal_planes",
      
      !is.na(Rrup_approx_np1_km) ~
        "nodal_plane_1_only",
      
      !is.na(Rrup_approx_np2_km) ~
        "nodal_plane_2_only",
      
      TRUE ~
        NA_character_
      
    ),
    
    rupture_centre_Rhypo_km = case_when(
      
      !is.na(rupture_centre_Repi_km) &
        !is.na(rupture_centre_depth_km) ~
        sqrt(
          rupture_centre_Repi_km^2 +
            rupture_centre_depth_km^2
        ),
      
      TRUE ~
        NA_real_
      
    )
    
  )


# ------------------------------------------------------------
# Select one convenience candidate for sensitivity analysis
# ------------------------------------------------------------

if (
  rrup_candidate_for_sensitivity ==
  "minimum"
) {
  
  iar_rrup$Rrup_approx_preferred_km <-
    iar_rrup$Rrup_approx_min_km
  
  
} else if (
  rrup_candidate_for_sensitivity ==
  "mean"
) {
  
  iar_rrup$Rrup_approx_preferred_km <-
    iar_rrup$Rrup_approx_mean_km
  
  
} else if (
  rrup_candidate_for_sensitivity ==
  "np1"
) {
  
  iar_rrup$Rrup_approx_preferred_km <-
    iar_rrup$Rrup_approx_np1_km
  
  
} else if (
  rrup_candidate_for_sensitivity ==
  "np2"
) {
  
  iar_rrup$Rrup_approx_preferred_km <-
    iar_rrup$Rrup_approx_np2_km
  
}


# ============================================================
# QUALITY AND AVAILABILITY FLAGS
# ============================================================

iar_rrup <- iar_rrup |>
  
  mutate(
    
    approximate_rrup_available =
      !is.na(
        Rrup_approx_preferred_km
      ),
    
    Rrup_distance_status = case_when(
      
      is.na(rupture_mw) ~
        "unavailable_missing_Mw",
      
      !rupture_scaling_eligible &
        rupture_scaling_extrapolated ~
        "unavailable_outside_scaling_range",
      
      !rupture_scaling_eligible ~
        "unavailable_scaling_ineligible",
      
      is.na(rupture_centre_longitude) |
        is.na(rupture_centre_latitude) |
        is.na(rupture_centre_depth_km) ~
        "unavailable_missing_rupture_centre",
      
      is.na(Rrup_approx_np1_km) &
        is.na(Rrup_approx_np2_km) ~
        "unavailable_missing_focal_geometry",
      
      !is.na(Rrup_approx_np1_km) &
        !is.na(Rrup_approx_np2_km) ~
        "available_both_nodal_planes",
      
      !is.na(Rrup_approx_np1_km) ~
        "available_nodal_plane_1_only",
      
      !is.na(Rrup_approx_np2_km) ~
        "available_nodal_plane_2_only",
      
      TRUE ~
        "unclassified"
      
    ),
    
    rrup_exceeds_centre_hypocentral =
      !is.na(Rrup_approx_min_km) &
      !is.na(rupture_centre_Rhypo_km) &
      Rrup_approx_min_km >
      rupture_centre_Rhypo_km +
      1e-6,
    
    distance_type =
      case_when(
        
        approximate_rrup_available ~
          paste0(
            "approximate_finite_rectangle_",
            rrup_candidate_for_sensitivity
          ),
        
        TRUE ~
          "hypocentral_only"
        
      )
    
  )


# ============================================================
# OPTIONALLY USE APPROXIMATE RRUP AS PRIMARY MODEL DISTANCE
# ============================================================

if (use_approx_rrup_as_primary_distance) {
  
  iar_rrup <- iar_rrup |>
    
    mutate(
      
      R_model_km = case_when(
        
        !is.na(
          Rrup_approx_preferred_km
        ) ~
          Rrup_approx_preferred_km,
        
        TRUE ~
          R_model_km
        
      ),
      
      R_model_type = case_when(
        
        !is.na(
          Rrup_approx_preferred_km
        ) ~
          paste0(
            "approximate_Rrup_",
            rrup_candidate_for_sensitivity
          ),
        
        TRUE ~
          R_model_type
        
      )
      
    )
  
}


# ============================================================
# DIAGNOSTICS
# ============================================================

rrup_availability_summary <- iar_rrup |>
  
  count(
    Rrup_distance_status,
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


rrup_event_summary <- iar_rrup |>
  
  group_by(
    publicid
  ) |>
  
  summarise(
    
    rupture_mw =
      first(rupture_mw),
    
    scaling_relation =
      first(scaling_relation),
    
    rupture_centre_source =
      first(rupture_centre_source),
    
    rupture_length_km =
      first(rupture_length_km),
    
    predicted_rupture_width_km =
      first(predicted_rupture_width_km),
    
    felt_cell_count =
      n(),
    
    cells_with_np1_distance =
      sum(
        !is.na(Rrup_approx_np1_km)
      ),
    
    cells_with_np2_distance =
      sum(
        !is.na(Rrup_approx_np2_km)
      ),
    
    cells_with_preferred_distance =
      sum(
        !is.na(Rrup_approx_preferred_km)
      ),
    
    minimum_Rrup_approx_km = if (
      all(
        is.na(
          Rrup_approx_preferred_km
        )
      )
    ) {
      NA_real_
    } else {
      min(
        Rrup_approx_preferred_km,
        na.rm = TRUE
      )
    },
    
    median_Rrup_approx_km = if (
      all(
        is.na(
          Rrup_approx_preferred_km
        )
      )
    ) {
      NA_real_
    } else {
      median(
        Rrup_approx_preferred_km,
        na.rm = TRUE
      )
    },
    
    maximum_Rrup_approx_km = if (
      all(
        is.na(
          Rrup_approx_preferred_km
        )
      )
    ) {
      NA_real_
    } else {
      max(
        Rrup_approx_preferred_km,
        na.rm = TRUE
      )
    },
    
    median_nodal_plane_range_km = if (
      all(
        is.na(
          Rrup_nodal_plane_range_km
        )
      )
    ) {
      NA_real_
    } else {
      median(
        Rrup_nodal_plane_range_km,
        na.rm = TRUE
      )
    },
    
    .groups = "drop"
    
  )


rrup_quality_flags <- iar_rrup |>
  
  filter(
    rrup_exceeds_centre_hypocentral |
      rupture_width_capped_np1 |
      rupture_width_capped_np2 |
      rupture_scaling_extrapolated
  ) |>
  
  select(
    publicid,
    cell_id = any_of("cell_id"),
    rupture_mw,
    scaling_relation,
    rupture_scaling_extrapolated,
    rupture_width_capped_np1,
    rupture_width_capped_np2,
    rupture_centre_Rhypo_km,
    Rrup_approx_np1_km,
    Rrup_approx_np2_km,
    Rrup_approx_min_km,
    rrup_exceeds_centre_hypocentral
  )


# ============================================================
# SAVE OUTPUTS
# ============================================================

arrow::write_parquet(
  iar_rrup,
  output_file
)


readr::write_csv(
  rrup_availability_summary,
  file.path(
    diagnostic_dir,
    "iar_rrup_availability_summary.csv"
  )
)


readr::write_csv(
  rrup_event_summary,
  file.path(
    diagnostic_dir,
    "iar_rrup_event_summary.csv"
  )
)


readr::write_csv(
  rrup_quality_flags,
  file.path(
    diagnostic_dir,
    "iar_rrup_quality_flags.csv"
  )
)


# ============================================================
# VERIFY OUTPUTS
# ============================================================

output_files <- c(
  rupture_geometry_file,
  output_file,
  file.path(
    diagnostic_dir,
    "iar_rrup_availability_summary.csv"
  ),
  file.path(
    diagnostic_dir,
    "iar_rrup_event_summary.csv"
  ),
  file.path(
    diagnostic_dir,
    "iar_rrup_quality_flags.csv"
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
      "One or more Script 6 outputs were not written correctly:\n",
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
  "\nApproximate rupture-distance processing complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Input observations:",
  format(
    nrow(iar),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Unique earthquakes:",
  format(
    n_distinct(
      iar_rrup$publicid,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Earthquakes with rupture dimensions:",
  format(
    sum(
      !is.na(
        rupture_events$rupture_length_km
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Observations with nodal-plane-1 Rrup:",
  format(
    sum(
      !is.na(
        iar_rrup$Rrup_approx_np1_km
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Observations with nodal-plane-2 Rrup:",
  format(
    sum(
      !is.na(
        iar_rrup$Rrup_approx_np2_km
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Observations with preferred approximate Rrup:",
  format(
    sum(
      !is.na(
        iar_rrup$Rrup_approx_preferred_km
      )
    ),
    big.mark = ","
  ),
  "\n"
)

cat(
  "Preferred candidate:",
  rrup_candidate_for_sensitivity,
  "\n"
)

cat(
  "Approximate Rrup used as primary model distance:",
  use_approx_rrup_as_primary_distance,
  "\n"
)

cat(
  "\nEvent-level rupture geometry:\n",
  normalizePath(
    rupture_geometry_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat(
  "\nFelt RAPID dataset with approximate Rrup:\n",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)


if (
  any(
    !is.na(
      iar_rrup$Rrup_approx_preferred_km
    )
  )
) {
  
  cat(
    "\nApproximate preferred Rrup summary:\n"
  )
  
  print(
    summary(
      iar_rrup$Rrup_approx_preferred_km
    )
  )
  
}