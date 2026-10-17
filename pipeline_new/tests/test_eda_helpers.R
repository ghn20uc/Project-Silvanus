#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(dirname(normalizePath(script))), "eda_helpers.R"))
suppressPackageStartupMessages(library(dplyr))
checks <- 0L
check <- function(value) { stopifnot(isTRUE(value)); checks <<- checks + 1L }
fails <- function(expression) inherits(tryCatch({force(expression); NULL}, error = identity), "error")
fixture <- data.frame(publicid = c("a", "a", "b"), cell_id = c("x", "y", "x"),
  origintime = as.POSIXct(c("2020-01-01", "2020-01-01", "2020-02-01"), tz = "UTC"),
  event_longitude = 172, event_latitude = -42, event_depth_km = 10,
  report_longitude = c(171, 172, 171), report_latitude = -43,
  reported_mmi = c(3, 3.5, 1), report_count = c(4, 6, 4), Mw = c(5, 5, 4),
  Repi_km = c(10, 20, 15), Rhypo_km = c(14, 22, 18), Vs30_m_s = c(200, NA, 300))
check(validate_cells(fixture, "fixture"))
# Current regression excludes MMI<3 and incomplete selected fields, but missing
# Rhypo is not itself a regression exclusion (it is a classification exclusion).
extra <- fixture[1, ]; extra$cell_id <- "z"; extra$Rhypo_km <- NA_real_
data <- rbind(fixture, extra)
minimum_reports_per_cell <- 4; minimum_reports_per_earthquake <- 4
prepared <- eval(regression_preparation)
check(identical(prepared$cell_id, c("x", "z")))
data$Repi_km[1] <- 0
check(identical(eval(regression_preparation)$cell_id, "z"))
source_text <- paste(readLines(file.path(dirname(dirname(normalizePath(script))),
  "06_models_modelmaps.Rmd"), warn = FALSE), collapse = "\n")
check(identical(notebook_assignment(source_text, "data_sd"), regression_preparation))
changed <- gsub("reported_mmi >= 3", "reported_mmi >= 2", source_text, fixed = TRUE)
check(!identical(notebook_assignment(changed, "data_sd"), regression_preparation))
changed <- sub("        cell_id,", "        Rhypo_km,", source_text, fixed = TRUE)
check(!identical(notebook_assignment(changed, "data_sd"), regression_preparation))
check(fails(validate_cells(rbind(fixture, fixture[1, ]), "duplicate")))
bad <- fixture; bad$Mw[2] <- 6
check(fails(validate_cells(bad, "conflicting event")))
bad <- fixture; bad$report_count[1] <- 4.5
check(fails(validate_cells(bad, "fractional count")))
bad <- fixture; bad$Vs30_m_s[1] <- Inf
check(fails(validate_cells(bad, "infinite value")))
bad <- fixture; bad$Repi_km <- NULL
check(fails(validate_cells(bad, "missing column")))
check(validate_split_ids(c("a", "b"), "c"))
check(fails(validate_split_ids("a", c("a", "b"))))
check(fails(validate_split_ids(c("a", "a"), "b")))
check(is.na(safe_spearman(1:9, 9:1)))
check(is.na(safe_spearman(rep(1, 10), 1:10)))
check(abs(safe_spearman(1:10, 10:1) + 1) < 1e-12)
check(numeric_summary(c(1, NA, 3))$median == 2)
check(numeric_summary(c(NA_real_, NA_real_))$n == 0)
# Publication must preserve unknown output files and replace only managed runs.
area <- tempfile("eda-test-"); dir.create(area)
dest <- file.path(area, "output"); dir.create(dest)
writeLines("keep", file.path(dest, "user.txt"))
stage <- file.path(area, "stage"); dir.create(stage)
check(fails(publish_eda(stage, dest)))
check(readLines(file.path(dest, "user.txt")) == "keep")
unlink(dest, recursive = TRUE)
writeLines("owned", file.path(stage, ".DATA309_EDA"))
writeLines("old", file.path(stage, "old.csv"))
write.csv(data.frame(file = c(".DATA309_EDA", "old.csv")), file.path(stage, "output_manifest.csv"), row.names = FALSE)
publish_eda(stage, dest)
next_stage <- file.path(area, "next"); dir.create(next_stage)
writeLines("owned", file.path(next_stage, ".DATA309_EDA"))
writeLines("new", file.path(next_stage, "new.csv"))
write.csv(data.frame(file = c(".DATA309_EDA", "new.csv")), file.path(next_stage, "output_manifest.csv"), row.names = FALSE)
publish_eda(next_stage, dest)
check(!file.exists(file.path(dest, "old.csv")) && file.exists(file.path(dest, "new.csv")))
unlink(area, recursive = TRUE)
cat(checks, "helper checks passed\n")
