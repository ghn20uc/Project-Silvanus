#!/usr/bin/env Rscript
# DATA309: exploratory analysis of the current final_data snapshot.
# Run from any directory: Rscript /path/to/Project-Silvanus/pipeline_new/DATA309_EDA.R
# See DATA309_EDA.md for inputs, sample definitions and output migration.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this entry point with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root <- dirname(script_dir)
source(file.path(script_dir, "eda_helpers.R"))
if ("--help" %in% args) {
  cat("Usage: Rscript DATA309_EDA.R [--data-dir PATH] [--output-dir PATH]\n",
      "Defaults: <repository>/final_data and <repository>/pipeline_new/eda_outputs\n")
  quit(status = 0)
}
options <- list("--data-dir" = file.path(root, "final_data"),
                "--output-dir" = file.path(script_dir, "eda_outputs"))
if (length(args) %% 2L) stop("Options require a path; use --help")
if (length(args)) for (i in seq(1, length(args), by = 2)) {
  if (!args[i] %in% names(options)) stop("Unknown option: ", args[i])
  options[[args[i]]] <- args[i + 1]
}
for (package in c("arrow", "dplyr", "tidyr", "ggplot2", "readr")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install the R package: ", package)
}
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
data_dir <- normalizePath(options[["--data-dir"]], mustWork = TRUE)
out_dir <- options[["--output-dir"]]
dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
out_dir <- file.path(normalizePath(dirname(out_dir)), basename(out_dir))
if (out_dir %in% c(root, script_dir, data_dir) || file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("Choose a dedicated EDA output directory")
}
inputs <- c(base = file.path(data_dir, "NZ_FR_IAR_dataset.parquet"),
  site = file.path(data_dir, "NZ_FR_IAR_final_dataset.parquet"),
  catalogue = file.path(data_dir, "geonet_catalogue_with_mw.parquet"),
  train_ids = file.path(script_dir, "data_raw/train_earthquake_ids.rds"),
  test_ids = file.path(script_dir, "data_raw/test_earthquake_ids.rds"),
  regression_source = file.path(script_dir, "06_models_modelmaps.Rmd"),
  classification_source = file.path(script_dir, "07_trees.qmd"),
  dataset_builder = file.path(root, "data_collection_pipeline_final/04_build_nz_fr_iar_dataset.R"),
  eda_source = file.path(script_dir, "DATA309_EDA.R"),
  eda_helpers = file.path(script_dir, "eda_helpers.R"))
if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
initial_hashes <- tools::md5sum(inputs)
base <- as.data.frame(arrow::read_parquet(inputs[["base"]]))
site <- as.data.frame(arrow::read_parquet(inputs[["site"]]))
catalogue <- as.data.frame(arrow::read_parquet(inputs[["catalogue"]]))
validate_cells(base, "Pre-site dataset")
validate_cells(site, "Final dataset")
require_columns(site, "Vs30_m_s", "Final dataset")
require_columns(catalogue, c("publicid", "origintime", "Mw"), "Catalogue")
if (anyNA(catalogue$publicid) || anyDuplicated(catalogue$publicid)) stop("Invalid catalogue IDs")
# A site enrichment must not change the underlying event-cell observations.
key <- function(x) paste(x$publicid, x$cell_id, sep = "\r")
idx <- match(key(site), key(base))
if (nrow(site) != nrow(base) || anyNA(idx) ||
    !isTRUE(all.equal(base[idx, names(base)], site[names(base)], check.attributes = FALSE))) {
  stop("Pre-site and final datasets disagree beyond site enrichment; inspect upstream data")
}
if (any(!site$publicid %in% catalogue$publicid)) stop("Final events missing from catalogue")
train_ids <- readRDS(inputs[["train_ids"]])
test_ids <- readRDS(inputs[["test_ids"]])
validate_split_ids(train_ids, test_ids)

# Fail visibly if the notebook constants this audit mirrors have changed.
reg_source <- paste(readLines(inputs[["regression_source"]], warn = FALSE), collapse = "\n")
class_source <- paste(readLines(inputs[["classification_source"]], warn = FALSE), collapse = "\n")
contracts <- c(
  grepl("minimum_reports_per_cell\\s*<-\\s*4\\b", reg_source, perl = TRUE),
  grepl("minimum_reports_per_earthquake\\s*<-\\s*4\\b", reg_source, perl = TRUE),
  grepl("set.seed\\(123\\)", reg_source),
  grepl("0\\.8\\s*\\*\\s*length", reg_source, perl = TRUE),
  grepl("filter\\(report_count >= 4\\)", class_source),
  grepl("levels = 3:8", class_source, fixed = TRUE))
if (!all(contracts)) stop("Model preparation constants changed; update the EDA sample rules before running")
if (!identical(notebook_assignment(reg_source, "data_sd"), regression_preparation)) {
  stop("Regression preparation changed; review the EDA fields and filtering rules before running")
}

# Match the current preparation rules without executing or editing model notebooks.
# Regression: >=4 reports per cell/event, complete selected fields, positive
# epicentral distance and Vs30, and MMI >=3. Match the approved notebook expression.
# The event threshold is redundant here because each retained cell already has >=4.
# Classification: complete cases in its selected columns, MMI>=3, ordered levels3:8,
# >=4 reports, then the saved IDs. A fractional MMI is not silently rounded.
class_fields <- c("publicid", "Mw", "Vs30_m_s", "Rhypo_km", "reported_mmi",
  "report_count", "event_depth_km", "Repi_km", "report_longitude", "report_latitude")
regression <- eval(regression_preparation, envir = list(data = site,
  minimum_reports_per_cell = 4, minimum_reports_per_earthquake = 4)) |>
  select(all_of(names(site)))
classification <- site[complete.cases(site[class_fields]), ] |>
  filter(report_count >= 4, reported_mmi >= 3)
if (any(!classification$reported_mmi %in% 3:8)) {
  stop("Current classification factor levels 3:8 cannot represent some MMI values; do not round")
}
if (!nrow(regression) || !nrow(classification)) stop("No eligible observations")
if (any(regression$Repi_km <= 0) || any(classification$Repi_km <= 0)) stop("Nonpositive log-distance input")
# Reconstruct the regression's split, preserving the input event order and R sample kind.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(123)
reg_ids <- unique(regression$publicid)
reg_train <- sample(reg_ids, floor(.8 * length(reg_ids)))
split_audit <- tibble(publicid = union(union(site$publicid, train_ids), test_ids)) |>
  mutate(in_final_data = publicid %in% site$publicid,
    regression_eligible = publicid %in% regression$publicid,
    classification_eligible = publicid %in% classification$publicid,
    saved_classification_split = case_when(publicid %in% train_ids ~ "Train",
      publicid %in% test_ids ~ "Test", TRUE ~ "Unassigned"),
    reconstructed_regression_split = case_when(publicid %in% reg_train ~ "Train",
      publicid %in% reg_ids ~ "Test", TRUE ~ "Ineligible"))
classification <- classification |> left_join(
  select(split_audit, publicid, saved_classification_split), by = "publicid")

# A temporary sibling directory prevents failed calculations leaving mixed runs.
stage <- tempfile(".eda-run-", tmpdir = dirname(out_dir))
dir.create(stage)
tryCatch({
writeLines("Owned outputs of pipeline_new/DATA309_EDA.R", file.path(stage, ".DATA309_EDA"))
write_table <- function(x, name) readr::write_csv(x, file.path(stage, paste0(name, ".csv")), na = "NA")
theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  plot.title.position = "plot", plot.caption = element_text(hjust = 0),
  plot.margin = margin(12, 18, 12, 12)))
save_plot <- function(plot, name, width = 9, height = 5.5) {
  ggsave(file.path(stage, paste0(name, ".png")), plot, width = width, height = height,
    dpi = 300, bg = "white")
}

stages <- list(pre_site = base, final_with_site = site,
  report_count_ge4 = filter(site, report_count >= 4), regression_eligible = regression,
  regression_train = filter(regression, publicid %in% reg_train),
  regression_test = filter(regression, !publicid %in% reg_train),
  classification_eligible = classification,
  classification_train = filter(classification, saved_classification_split == "Train"),
  classification_test = filter(classification, saved_classification_split == "Test"),
  classification_unassigned = filter(classification, saved_classification_split == "Unassigned"))
overview <- bind_rows(lapply(names(stages), function(nm) {
  x <- stages[[nm]]
  tibble(stage = nm, earthquake_cells = nrow(x), earthquakes = n_distinct(x$publicid),
    distinct_cell_ids = n_distinct(x$cell_id),
    reporting_coordinates = nrow(unique(x[c("report_longitude", "report_latitude")])),
    report_contributions = sum(x$report_count),
    earliest_utc = if (nrow(x)) format(min(x$origintime), tz = "UTC", usetz = TRUE) else NA_character_,
    latest_utc = if (nrow(x)) format(max(x$origintime), tz = "UTC", usetz = TRUE) else NA_character_)
}))
write_table(overview, "00_dataset_overview")
write_table(tibble(metric = c("catalogue_earthquakes", "catalogue_with_Mw", "final_earthquakes"),
  value = c(nrow(catalogue), sum(!is.na(catalogue$Mw)), n_distinct(site$publicid))), "00_catalogue_overview")
write_table(split_audit, "01_earthquake_split_audit")
write_table(count(classification, saved_classification_split, reported_mmi, name = "earthquake_cells"),
  "01_split_mmi_counts")
write_table(regression |> mutate(split = if_else(publicid %in% reg_train, "Train", "Test")) |>
  count(split, reported_mmi, name = "earthquake_cells"), "01_regression_split_mmi_counts")
write_table(count(filter(split_audit, classification_eligible), saved_classification_split,
  reconstructed_regression_split, name = "earthquakes"), "01_split_comparison")

mmi <- bind_rows(lapply(names(stages), function(nm) {
  count(stages[[nm]], reported_mmi, name = "earthquake_cells") |>
    mutate(stage = nm, percent = 100 * earthquake_cells / sum(earthquake_cells))
}))
write_table(mmi, "02_mmi_distribution")
plot_mmi <- filter(mmi, stage %in% c("final_with_site", "regression_eligible", "classification_eligible"))
save_plot(ggplot(plot_mmi, aes(factor(reported_mmi), earthquake_cells)) +
  geom_col(fill = "#217c83") + facet_wrap(~stage, ncol = 1, scales = "free_y") +
  labs(title = "Reported MMI across sample definitions", x = "Median reported MMI",
    y = "Earthquake-cell observations", caption = "Both model preparations require MMI >= 3. The final dataset retains the lower values."),
  "02_mmi_distribution", height = 8)

missingness <- bind_rows(lapply(c("pre_site", "final_with_site", "report_count_ge4"), function(nm) {
  x <- stages[[nm]]
  tibble(stage = nm, variable = names(x), missing = vapply(x, function(v) sum(is.na(v)), integer(1)),
    rows = nrow(x)) |> mutate(percent_missing = 100 * missing / rows)
}))
write_table(missingness, "03_missingness")
thresholds <- bind_rows(lapply(c(1, 2, 3, 4, 5, 10, 20, 50, 100), function(threshold) {
  x <- site |> filter(report_count >= threshold)
  tibble(minimum_reports = threshold, earthquake_cells = nrow(x), earthquakes = n_distinct(x$publicid),
    percent_all_cells = 100 * nrow(x) / nrow(site), complete_cells = sum(complete.cases(x)),
    cells_mmi_ge5 = sum(x$reported_mmi >= 5, na.rm = TRUE),
    cells_mmi_ge7 = sum(x$reported_mmi >= 7, na.rm = TRUE))
}))
write_table(thresholds, "04_report_threshold_retention")
save_plot(ggplot(thresholds, aes(minimum_reports, percent_all_cells)) + geom_line(colour = "#217c83") +
  geom_point() + scale_x_log10(breaks = thresholds$minimum_reports) +
  geom_vline(xintercept = 4, linetype = 2, colour = "#ad503b") +
  labs(title = "Observations retained by the report-count threshold", x = "Minimum reports per earthquake-cell",
    y = "% of final dataset retained", caption = "Before complete-case filtering. The dashed line marks the current threshold, not an estimated optimum."),
  "04_report_threshold_retention")

events <- distinct(site, publicid, origintime, Mw, event_depth_km, event_longitude, event_latitude)
event_totals <- site |> group_by(publicid) |> summarise(cells = n(),
  report_contributions = sum(report_count), .groups = "drop")
write_table(left_join(events, event_totals, by = "publicid"), "05_event_coverage")
write_table(events |> mutate(year = as.integer(format(origintime, "%Y", tz = "UTC"))) |>
  count(year, name = "earthquakes"), "05_events_by_year")
# Earthquake attributes are counted once per event; the other summaries use cells.
variables <- list(Mw = events$Mw, event_depth_km = events$event_depth_km,
  report_count = site$report_count, cells_per_earthquake = event_totals$cells,
  reports_per_earthquake = event_totals$report_contributions,
  Repi_km = regression$Repi_km, Rhypo_km = regression$Rhypo_km, Vs30_m_s = regression$Vs30_m_s)
units <- c("earthquake (final)", "earthquake (final)", "earthquake-cell (final)",
  "earthquake (final)", "earthquake (final)", rep("earthquake-cell (regression eligible)", 3))
write_table(bind_rows(lapply(seq_along(variables), function(i) {
  cbind(variable = names(variables)[i], unit = units[i], numeric_summary(variables[[i]]))
})), "06_numeric_summary")
for (field in c("Mw", "event_depth_km")) {
  save_plot(ggplot(events, aes(x = .data[[field]])) + geom_histogram(bins = 40, fill = "#217c83", colour = "white") +
    labs(title = paste("Earthquake-level distribution:", field), x = field, y = "Earthquakes",
      caption = "One observation per earthquake represented in the final dataset."), paste0("06_", field, "_distribution"))
}
save_plot(ggplot(site, aes(report_count)) + geom_histogram(bins = 45, fill = "#217c83", colour = "white") +
  scale_x_log10() + labs(title = "Public reports per earthquake-cell", x = "Report count (log scale)",
    y = "Earthquake-cell observations", caption = "Report contributions are not a count of unique people across earthquakes."), "06_report_count_distribution")
save_plot(ggplot(regression, aes(Vs30_m_s)) + geom_histogram(bins = 40, fill = "#217c83", colour = "white") +
  labs(title = "Site conditions in the regression-eligible sample", x = "Vs30 (m/s)", y = "Earthquake-cell observations",
    caption = "Locations can recur across earthquakes; this is not a distribution of unique sites."), "06_Vs30_distribution")

# Descriptive associations: no causal interpretation or independent-cell p-values.
relations <- regression |> mutate(magnitude_band = cut(Mw, c(-Inf, 4, 5, 6, Inf),
  right = FALSE, labels = c("Mw < 4", "4 <= Mw < 5", "5 <= Mw < 6", "Mw >= 6")),
  distance_band = cut(Repi_km, c(0, 10, 20, 50, 100, 200, 500, Inf), include.lowest = TRUE))
attenuation <- relations |> group_by(magnitude_band, distance_band) |> summarise(
  earthquake_cells = n(), earthquakes = n_distinct(publicid), distance_median_km = median(Repi_km),
  mmi_q25 = quantile(reported_mmi, .25), mmi_median = median(reported_mmi),
  mmi_q75 = quantile(reported_mmi, .75), .groups = "drop")
write_table(attenuation, "07_distance_by_magnitude")
save_plot(ggplot(attenuation, aes(distance_median_km, mmi_median)) +
  geom_linerange(aes(ymin = mmi_q25, ymax = mmi_q75), colour = "#76aeb0", linewidth = 1) +
  geom_line(colour = "#217c83") + geom_point(aes(size = earthquake_cells)) +
  facet_wrap(~magnitude_band) + scale_x_log10() +
  labs(title = "Shaking and epicentral distance, grouped by magnitude", x = "Median distance within bin (km, log scale)",
    y = "Median reported MMI", size = "Cells", caption = "Regression-eligible cells. Bars show the interquartile range, not uncertainty intervals.\nDifferent earthquakes contribute at different distances; these are descriptive summaries."),
  "07_distance_by_magnitude", height = 7)
within_event <- regression |> group_by(publicid) |> summarise(earthquake_cells = n(),
  Mw = first(Mw), rho_Repi = safe_spearman(Repi_km, reported_mmi),
  rho_Rhypo = safe_spearman(Rhypo_km, reported_mmi), rho_Vs30 = safe_spearman(Vs30_m_s, reported_mmi),
  .groups = "drop")
write_table(within_event, "08_within_event_associations")
write_table(tibble(variable = c("Repi_km", "Rhypo_km", "Mw", "Vs30_m_s"),
  pooled_spearman = vapply(c("Repi_km", "Rhypo_km", "Mw", "Vs30_m_s"), function(v) {
    safe_spearman(regression[[v]], regression$reported_mmi)
  }, numeric(1))), "08_pooled_associations")
vs_bins <- regression |> mutate(Vs30_band = cut(Vs30_m_s, c(0, 180, 360, 760, Inf), right = FALSE))
write_table(vs_bins |> group_by(Vs30_band) |> summarise(earthquake_cells = n(),
  earthquakes = n_distinct(publicid), mmi_median = median(reported_mmi), .groups = "drop"), "09_site_conditions")
save_plot(ggplot(vs_bins, aes(Vs30_band, reported_mmi)) + geom_boxplot(fill = "#b3d7d7", outlier.alpha = .05) +
  labs(title = "Reported intensity across Vs30 ranges", x = "Vs30 range (m/s)", y = "Median reported MMI",
    caption = "Regression-eligible cells. Magnitude, distance and event composition also vary between groups."), "09_site_conditions")
# Select examples by coverage, not by strength of correlation.
example_ids <- within_event |> arrange(desc(earthquake_cells), publicid) |> head(4) |> pull(publicid)
examples <- regression |> filter(publicid %in% example_ids) |>
  mutate(event = paste0(publicid, " (Mw ", Mw, ")"))
save_plot(ggplot(examples, aes(Repi_km, reported_mmi)) + geom_point(alpha = .10, size = .6, colour = "#217c83") +
  scale_x_log10() + facet_wrap(~event) + labs(title = "Distance and intensity within individual earthquakes",
    x = "Epicentral distance (km, log scale)", y = "Median reported MMI",
    caption = "Four earthquakes with the most regression-eligible cells; all their eligible cells are shown."),
  "10_event_examples", height = 7)

# Existing local outline only: no API calls or GIS downloads.
outline_file <- file.path(script_dir, "nzmap/regional-council-2025-clipped.gpkg")
map_status <- "Skipped: local outline or optional sf package unavailable"
if (file.exists(outline_file) && requireNamespace("sf", quietly = TRUE)) {
  outline <- sf::st_read(outline_file, quiet = TRUE) |> sf::st_transform(4326)
  locations <- regression |> count(report_longitude, report_latitude, name = "earthquake_cells")
  map <- ggplot() + geom_sf(data = outline, fill = "#f0f0ed", colour = "#bdbdb8", linewidth = .15) +
    geom_point(data = locations, aes(report_longitude, report_latitude, colour = log10(earthquake_cells)), size = .5) +
    scale_colour_viridis_c(name = "log10(cells)") +
    coord_sf(xlim = c(165, 180), ylim = c(-48, -33), expand = FALSE) +
    labs(title = "Coverage of reporting locations", x = "Longitude", y = "Latitude",
      caption = "Regression-eligible earthquake-cells at each coordinate.\nMain New Zealand extent; all coordinates remain in the accompanying table.")
  save_plot(map, "11_reporting_coverage", width = 7, height = 8)
  write_table(locations, "11_reporting_locations")
  inputs <- c(inputs, map_outline = outline_file)
  map_status <- "Generated from the repository's local regional-council outline"
}

availability <- tibble(analysis = c("Within-cell mean, mode, spread and skewness",
  "Median precision or mean/median comparison", "Fault-style comparisons", "Geographic coverage"),
  status = c(rep("Unavailable", 3), if (startsWith(map_status, "Generated")) "Generated" else "Skipped"),
  reason = c("Current files contain only cell median and report count, not individual MMI frequencies.",
    "Report count alone does not identify within-cell variability or sampling uncertainty.",
    "The current final dataset has no fault_style field.", map_status))
write_table(availability, "12_analysis_availability")

unassigned <- sum(classification$saved_classification_split == "Unassigned")
comp <- filter(split_audit, classification_eligible)
disagree <- sum(comp$saved_classification_split != comp$reconstructed_regression_split)
same_eligible <- setequal(key(regression), key(classification))
reg_old_path <- grepl('"data_processed"', reg_source, fixed = TRUE)
class_current_path <- grepl('"../final_data/NZ_FR_IAR_final_dataset.parquet"', class_source, fixed = TRUE)
builder_source <- paste(readLines(inputs[["dataset_builder"]], warn = FALSE), collapse = "\n")
strict_mw35 <- grepl("minimum_moment_magnitude\\s*<-\\s*3\\.5\\b", builder_source, perl = TRUE) &&
  grepl("Mw\\s*>\\s*minimum_moment_magnitude", builder_source, perl = TRUE)
mw_conflicts <- if (strict_mw35) sum(site$Mw <= 3.5, na.rm = TRUE) else NA_integer_
write_table(tibble(check = c("regression_input_location", "classification_input_location",
    "saved_split_coverage", "split_agreement", "sample_definition", "upstream_magnitude_filter"),
  status = c(if (reg_old_path) "Review" else "Check source",
    if (class_current_path) "Pass" else "Check source",
    if (unassigned) "Review" else "Pass", if (disagree) "Review" else "Pass", "Documented",
    if (is.na(mw_conflicts)) "Check source" else if (mw_conflicts) "Review" else "Pass"),
  detail = c(if (reg_old_path) "Regression source still reads data_processed; EDA reads final_data." else "Review the current regression input expression.",
    if (class_current_path) "Classification reads ../final_data from pipeline_new." else "Review the current classification input expression.",
    sprintf("%s eligible classification cells have no saved training/test earthquake ID.", unassigned),
    sprintf("%s classification-eligible earthquakes differ between the saved split and seed-123 regression reconstruction.", disagree),
    sprintf("Regression: %s cells. Classification before saved-ID selection: %s cells. Identical eligible event-cell keys: %s; partitions are audited separately.", nrow(regression), nrow(classification), same_eligible),
    if (is.na(mw_conflicts)) "Builder threshold changed; inspect its current magnitude rule." else
      sprintf("Builder requires Mw > 3.5, but current final data contain %s rows at or below 3.5. EDA retains the supplied data and does not resolve this policy difference.", mw_conflicts))),
  "13_pipeline_audit")

# Include source hashes and actual data dates; do not reuse old report counts/dates.
if (!identical(initial_hashes, tools::md5sum(names(initial_hashes)))) stop("Inputs changed during EDA")
manifest <- tibble(role = names(inputs), file = vapply(inputs, function(p) {
  if (startsWith(p, paste0(root, "/"))) substring(p, nchar(root) + 2L) else p
}, character(1)), bytes = file.info(inputs)$size, md5 = unname(tools::md5sum(inputs)))
write_table(manifest, "input_manifest")
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(stage, "sessionInfo.txt"))
run_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)
writeLines(c("# DATA309 EDA outputs", "", paste("Generated:", run_time), "",
  "Source: current final_data Parquet files. Observation: one earthquake-cell pair; response: median reported MMI.",
  "", "## Current sample", "",
  sprintf("Final data: %s earthquake-cells across %s earthquakes, from %s to %s (UTC).",
    nrow(site), n_distinct(site$publicid), min(site$origintime), max(site$origintime)),
  sprintf("Regression preparation retains %s earthquake-cells across %s earthquakes with at least four reports per cell, complete selected fields, positive Repi/Vs30 and MMI >=3.",
    nrow(regression), n_distinct(regression$publicid)),
  sprintf("Classification requires complete selected fields, at least four reports and MMI levels 3-8: %s eligible cells; %s are not covered by the saved IDs.", nrow(classification), unassigned),
  sprintf("Identical eligible event-cell keys before partitioning: %s. This does not establish identical train/test sets. See 01_split_comparison.csv and 13_pipeline_audit.csv.", same_eligible),
  "The current final regression predicts raw MMI with linear Vs30; this EDA describes observations and does not reproduce model performance.",
  "", "## Reading the results", "",
  "00: sample flow and catalogue coverage. 01: saved-ID coverage and reconstructed split comparison.",
  "02-06: MMI, missingness, threshold retention, event coverage and distributions.",
  "07-10: descriptive distance/magnitude/site associations and within-event examples.",
  "11: optional local coverage map. 12: unavailable analyses. 13: model integration issues outside EDA scope.",
  "All percentages in 02 use their own stage as denominator; 04 uses all final earthquake-cells.",
  "Vs30 and response summaries count earthquake-cells, so repeated locations/events receive repeated weight.",
  "Cell IDs in the current inputs include the earthquake identifier; distinct_cell_ids is not a count of unique physical locations. reporting_coordinates counts distinct coordinate pairs.",
  "Spearman correlations require at least ten observations and variation in both variables. NA means not estimable, not zero.",
  "Full-sample EDA is descriptive; it is not held-out performance validation or evidence for tuning against test observations.",
  "", "## Inputs and reproducibility", "",
  "input_manifest.csv identifies input and source-code hashes. output_manifest.csv identifies this run's products.",
  "No models are fitted, no data are downloaded, and no input datasets or saved IDs are overwritten.",
  "13_pipeline_audit.csv checks model input references, saved partitions and the upstream Mw threshold against this snapshot.",
  "See ../DATA309_EDA.md for commands, dependencies, migration and limitations."), file.path(stage, "README.md"))
share34 <- 100 * mean(regression$reported_mmi %in% c(3, 4))
writeLines(c("# Data-derived EDA text", "", "Regenerate this text when inputs change; check the sample used by your models before copying it into the report.", "",
  sprintf("The final dataset contains %s earthquake-cell observations. Requiring at least four reports per cell, complete selected fields, positive epicentral distance and Vs30, and MMI >=3 retains %s observations from %s earthquakes for regression.", nrow(site), nrow(regression), n_distinct(regression$publicid)),
  sprintf("In this sample, MMI 3 and 4 together account for %.2f%% of observations, and %s observations have MMI 5 or higher. These are cell-level medians, not individual report ratings.", share34, sum(regression$reported_mmi >= 5)),
  sprintf("The final data contain %s missing Vs30 values. Classification uses MMI 3-8 and saved event partitions; %s otherwise eligible observations are outside those partitions. Identical eligible event-cell keys before splitting: %s.", sum(is.na(site$Vs30_m_s)), unassigned, same_eligible),
  sprintf("The reconstructed regression split contains %s training observations across %s earthquakes and %s test observations across %s earthquakes. These partitions differ from the saved classification partitions.",
    sum(regression$publicid %in% reg_train), length(reg_train),
    sum(!regression$publicid %in% reg_train), length(setdiff(reg_ids, reg_train))),
  "The report-count retention analysis describes the sample-size trade-off; it does not establish that four reports is an optimal precision threshold.",
  "Individual MMI frequencies are not included in the current inputs, so within-cell dispersion and comparisons of mean, median and mode cannot be estimated from these files."), file.path(stage, "report_text.md"))
files <- list.files(stage, all.files = TRUE, no.. = TRUE)
write_table(tibble(file = files, bytes = file.info(file.path(stage, files))$size,
  md5 = unname(tools::md5sum(file.path(stage, files)))), "output_manifest")
publish_eda(stage, out_dir)
cat("EDA complete:", out_dir, "\n")
print(overview)
cat("Classification cells without saved IDs:", unassigned, "\n")
}, finally = if (dir.exists(stage)) unlink(stage, recursive = TRUE))
