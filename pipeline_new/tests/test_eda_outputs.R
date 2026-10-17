#!/usr/bin/env Rscript
# Independently recompute sample totals from raw inputs using base R.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
pipeline <- dirname(dirname(normalizePath(script)))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else file.path(pipeline, "eda_outputs")
data_dir <- if (length(args) >= 2) args[2] else file.path(dirname(pipeline), "final_data")
read <- function(name) read.csv(file.path(out, paste0(name, ".csv")), stringsAsFactors = FALSE)
site <- as.data.frame(arrow::read_parquet(file.path(data_dir, "NZ_FR_IAR_final_dataset.parquet")))
reg_fields <- c("publicid", "cell_id", "reported_mmi", "Mw", "event_depth_km", "Repi_km", "Vs30_m_s")
reg <- site[which(site$report_count >= 4 & complete.cases(site[reg_fields]) &
  site$Repi_km > 0 & site$Vs30_m_s > 0 & site$reported_mmi >= 3), ]
fields <- c("publicid", "Mw", "Vs30_m_s", "Rhypo_km", "reported_mmi", "report_count",
  "event_depth_km", "Repi_km", "report_longitude", "report_latitude")
cls <- site[which(site$report_count >= 4 & site$reported_mmi >= 3 & complete.cases(site[fields])), ]
train <- readRDS(file.path(pipeline, "data_raw/train_earthquake_ids.rds"))
test <- readRDS(file.path(pipeline, "data_raw/test_earthquake_ids.rds"))
expected <- c(final_with_site = nrow(site), regression_eligible = nrow(reg), classification_eligible = nrow(cls),
  classification_train = sum(cls$publicid %in% train), classification_test = sum(cls$publicid %in% test),
  classification_unassigned = sum(!cls$publicid %in% c(train, test)))
overview <- read("00_dataset_overview")
stopifnot(identical(as.numeric(expected), as.numeric(overview$earthquake_cells[match(names(expected), overview$stage)])))
stopifnot(sum(expected[c("classification_train", "classification_test", "classification_unassigned")]) == nrow(cls))
stopifnot(overview$reporting_coordinates[overview$stage == "final_with_site"] ==
  nrow(unique(site[c("report_longitude", "report_latitude")])))
mmi <- subset(read("02_mmi_distribution"), stage == "regression_eligible")
actual <- table(reg$reported_mmi)
stopifnot(identical(as.numeric(actual), as.numeric(mmi$earthquake_cells[match(names(actual), mmi$reported_mmi)])),
  abs(sum(mmi$percent) - 100) < 1e-8)
thresholds <- read("04_report_threshold_retention")
for (i in seq_len(nrow(thresholds))) {
  idx <- which(site$report_count >= thresholds$minimum_reports[i])
  stopifnot(thresholds$earthquake_cells[i] == length(idx),
    thresholds$cells_mmi_ge7[i] == sum(site$reported_mmi[idx] >= 7))
}
missing <- subset(read("03_missingness"), stage == "final_with_site")
stopifnot(all(missing$missing == colSums(is.na(site))[missing$variable]))
audit <- read("01_earthquake_split_audit")
set.seed(123)
ids <- unique(reg$publicid)
rt <- sample(ids, floor(.8 * length(ids)))
stopifnot(overview$earthquake_cells[overview$stage == "regression_train"] == sum(reg$publicid %in% rt),
  overview$earthquake_cells[overview$stage == "regression_test"] == sum(!reg$publicid %in% rt))
reg_counts <- read("01_regression_split_mmi_counts")
for (split in c("Train", "Test")) {
  selected <- if (split == "Train") reg$publicid %in% rt else !reg$publicid %in% rt
  counts <- table(reg$reported_mmi[selected])
  reported <- reg_counts[reg_counts$split == split, ]
  stopifnot(all(as.numeric(counts) == reported$earthquake_cells[match(names(counts), reported$reported_mmi)]))
}
stopifnot(setequal(audit$publicid[audit$reconstructed_regression_split == "Train"], rt),
  setequal(audit$publicid[audit$saved_classification_split == "Test"], test))
event_counts <- table(reg$publicid)
within <- read("08_within_event_associations")
for (id in head(names(sort(event_counts, decreasing = TRUE)), 4)) {
  x <- reg[reg$publicid == id, ]
  stopifnot(abs(within$rho_Repi[within$publicid == id] - cor(x$Repi_km, x$reported_mmi, method = "spearman")) < 1e-12)
}
manifest <- read("output_manifest")
inputs <- read("input_manifest")
input_paths <- ifelse(grepl("^/|^[A-Za-z]:", inputs$file), inputs$file, file.path(dirname(pipeline), inputs$file))
stopifnot(identical(unname(tools::md5sum(input_paths)), inputs$md5))
paths <- file.path(out, manifest$file)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size == manifest$bytes),
  identical(unname(tools::md5sum(paths)), manifest$md5),
  setequal(list.files(out, all.files = TRUE, no.. = TRUE), c(manifest$file, "output_manifest.csv")))
pngs <- paths[grepl("\\.png$", paths)]
stopifnot(length(pngs) >= 10)
for (path in pngs) {
  bytes <- readBin(path, "raw", n = 24)
  stopifnot(identical(bytes[1:8], as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))))
  dimensions <- sapply(c(17, 21), function(i) sum(as.integer(bytes[i:(i+3)]) * 256^(3:0)))
  stopifnot(all(dimensions >= 1500))
}
cat("Output validation passed: sample totals, MMI counts, thresholds, missingness, splits, correlations, manifests and", length(pngs), "PNG files\n")
