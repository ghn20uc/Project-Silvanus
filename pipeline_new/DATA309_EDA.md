# DATA309 EDA

`DATA309_EDA.R` describes the current repository data and audits the samples used by
`06_models_modelmaps.Rmd` and `07_trees.qmd`. It does not fit models or change either
notebook, the input data, or saved training/test IDs.

## Run

With R 4.5 and packages `arrow`, `dplyr` (>=1.1), `tidyr`, `ggplot2` and `readr`:

```sh
Rscript pipeline_new/DATA309_EDA.R
Rscript pipeline_new/tests/test_eda_helpers.R
Rscript pipeline_new/tests/test_eda_outputs.R
```

Install `sf` to also generate the geographic coverage map from the existing
`pipeline_new/nzmap/regional-council-2025-clipped.gpkg`. No network requests are
made by EDA. The script can be invoked by absolute path from any working directory.
No `setwd()`, IDE project state, `here()` configuration or preloaded workspace is needed.
Use the Rscript entry point rather than selecting all lines and clicking Run.
In the Windows RStudio Console, with the repository as the working directory:

```r
system2(file.path(R.home("bin"), "Rscript.exe"),
        args = shQuote("pipeline_new/DATA309_EDA.R"))
```

```sh
Rscript /path/to/Project-Silvanus/pipeline_new/DATA309_EDA.R \
  --data-dir /path/to/Project-Silvanus/final_data \
  --output-dir /tmp/data309-eda
Rscript /path/to/Project-Silvanus/pipeline_new/tests/test_eda_outputs.R \
  /tmp/data309-eda /path/to/Project-Silvanus/final_data
```

Default output: `pipeline_new/eda_outputs/`. The output directory is dedicated to
this script: a successful run replaces its previous managed products, including
optional products that are no longer available. Existing unknown directories or
extra user files are refused rather than overwritten. Calculations finish in a
temporary sibling directory before publication, preserving the previous successful
run on calculation failure. Do not manually edit generated CSVs or figures.

## Current inputs and migration

| Previous assumption | Current source / treatment |
| --- | --- |
| `data_processed/NZ_FR_IAR_with_site_conditions.parquet` | `final_data/NZ_FR_IAR_final_dataset.parquet` |
| Pre-site observations under `data_raw` | `final_data/NZ_FR_IAR_dataset.parquet`; validated against the final event-cell records |
| `catalogue_with_felt_mt_strong.parquet` | `final_data/geonet_catalogue_with_mw.parquet` |
| `mt_Mw`, `has_felt_rapid`, `vs30_status`, `fault_style` | Use `Mw` and observed `Vs30_m_s` missingness; do not infer absent fields |
| Individual reports / `felt_rapid_mmi_counts.parquet` | Not present in the current inputs; within-cell moments and median precision are unavailable |
| Working-directory-relative `eda_outputs/` | Stable path `pipeline_new/eda_outputs/` |
| Old fixed narrative and historical counts | Generated `eda_outputs/report_text.md`, with current sample definitions |

The original numbered output names have been replaced with the output index below.
A repository search found no other scripts reading the old EDA output filenames.
Old manually copied figures should be replaced using this index; legacy files in
other directories are not deleted. In particular, the old `06*` mean/median/mode
and precision plots, `07*` within-cell skewness plots and `10b_vs30_status_summary`
are retired, not regenerated using unsupported assumptions.

| Current prefix | Contents |
| --- | --- |
| `00` | Dataset stages, time ranges, catalogue coverage |
| `01` | Saved earthquake-ID coverage, class counts, comparison with reconstructed regression split |
| `02` | MMI distribution for each explicitly named sample |
| `03` | Missingness before/after report-count filtering |
| `04` | Report-count threshold retention, including high-MMI cells |
| `05` | Event attributes, report/cell coverage and annual event counts |
| `06` | Numeric summaries and distributions (Mw, depth, report count, Vs30) |
| `07` | Distance summaries within magnitude bands |
| `08` | Pooled and within-event Spearman associations |
| `09` | Site-condition summaries and descriptive boxplots |
| `10` | Four events with the most eligible cells, selected without using correlation strength |
| `11` | Geographic coverage; optional local map and coordinate table |
| `12` | Availability of analyses that require additional inputs |
| `13` | Integration issues in model notebooks, reported without modifying them |

Figures are 300-dpi PNGs. `input_manifest.csv` records source and data MD5 hashes
(for change detection), `output_manifest.csv` records every generated product except
itself, and `sessionInfo.txt` records R/package versions. `README.md` and
`report_text.md` are generated from the actual input snapshot. Generation time is
UTC; the observed date range comes from `origintime`, not filesystem timestamps.

## Samples and interpretation

- **Final dataset:** all event-cell rows in the enriched Parquet file, including
  missing Vs30 values. A row is an earthquake-cell pair, not an individual person.
- **Regression eligible:** at least four reports per cell, at least four per event,
  then complete `publicid`, `cell_id`, `reported_mmi`, `Mw`, `event_depth_km`,
  `Repi_km` and `Vs30_m_s`; require positive `Repi_km`/`Vs30_m_s` and MMI >=3.
  The event threshold is redundant after this cell threshold. This matches the
  updated regression notebook; lower MMI values remain in the full-data overview.
  Missingness in fields outside this list is not an extra regression exclusion.
- **Classification eligible:** complete cases in the notebook's selected fields,
  at least four reports, MMI >=3 and representable by its ordered levels 3-8.
  Unexpected fractional/out-of-range values cause an error; they are not rounded.
- **Classification train/test:** intersect eligible rows with the existing RDS
  event IDs. Eligible events outside both lists are explicitly **Unassigned**.
- **Reconstructed regression split:** R's seed 123, floor(0.8 * number of eligible
  events), preserving event order in the input. This audits current preparation
  logic; it does not establish which split produced previously committed models.
  The overview and `01_regression_split_mmi_counts.csv` report this split separately
  from the classification counts. Eligible event-cell keys currently match across
  the two methods, but their actual training/test partitions do not.

The current final regression uses raw MMI and untransformed Vs30, with both linear
and log10 epicentral distance. EDA charts describe observed MMI; they do not apply
the earlier model's Box-Cox transformation or claim current prediction accuracy.

The EDA uses the full data for descriptive coverage and associations. It is not
independent validation of a trained model, and full-sample trends should not be
used to tune against a held-out test set. Earthquake-level distributions count each
earthquake once; cell-level distributions weight earthquakes/locations by their
number of observations. The map counts event-cell observations at each coordinate.
Current `cell_id` values include the earthquake identifier, so they must not be
counted as unique physical sites across events. The overview separately reports
distinct coordinate pairs; matching coordinates do not establish identical cell
boundaries across events.
Report-count sums represent report contributions, not unique people.

Distance-bin bars show the IQR of cell medians, not a confidence interval or
within-cell variation. Vs30 relationships and pooled correlations do not isolate
site effects from magnitude, distance or earthquake composition. Within-event
Spearman summaries require ten cells and variation in both variables. `NA` means
not estimable, not no association.

Median plus report count cannot recover the original mean, mode, standard
deviation, skewness or uncertainty. More reports do not by themselves demonstrate
agreement or establish an optimal threshold of four. Those analyses require a
separately validated individual-report or per-MMI-frequency dataset.

## Model integration checks

The regression notebook still references `data_processed`; classification now
reads `../final_data` from `pipeline_new`. EDA reads `final_data` directly.
It audits saved partitions and the upstream magnitude boundary, records input/source hashes,
and reports discrepancies; it does not rewrite or execute the model notebooks.
Sample and split counts describe what their preparation rules imply for **current**
data, not previously fitted results. Review `13_pipeline_audit.csv` before claiming
identical training/test samples or reusing old model metrics with new data.

When model preparation rules change, review this script's matching rules and tests.
The parsed `data_sd` preparation expression must match the approved expression
in `eda_helpers.R`; changing the selected fields, MMI cutoff or filter order fails
before output publication. Threshold, factor-level and seed guards also apply.
Schema, event-cell uniqueness,
finite values, event attributes, enrichment consistency and split overlap are
validated before outputs are published.

The upstream builder currently requires Mw >3.5 while the supplied data include
Mw=3.5. This is recorded as a review item; EDA does not remove these observations
or silently decide a new threshold. Integration questions for the team are listed
in [DATA309_EDA_integration.md](DATA309_EDA_integration.md).
