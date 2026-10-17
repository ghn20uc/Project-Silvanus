# EDA integration questions

Reviewed against repository revision `b2f44ca` on 2026-10-18 (New Zealand date).
These are questions for the relevant pipeline/model owners, not changes made by
EDA. For regenerated counts and source hashes, use `eda_outputs/13_pipeline_audit.csv`
and `eda_outputs/input_manifest.csv`. This note records the inspected snapshot;
it is not automatically regenerated when those other files change.

## Data generation

1. **Which magnitude boundary is intended?**
   `data_collection_pipeline_final/04_build_nz_fr_iar_dataset.R` requires Mw >3.5,
   but the supplied data include 11,420 cells from 87 earthquakes with Mw exactly
   3.5. EDA retains the supplied data; changing to >=3.5 requires an explicit
   project decision and matching code/data, not merely changing documentation.

2. **Which scripts generated the latest Parquet files?**
   Comparing upload `9104d03` with `56a8760`: the catalogue is byte-identical;
   the base dataset gains 11,420 Mw=3.5 rows and loses 459 other rows. All retained
   base rows have unchanged field values. Among shared final rows, 14,937 Vs30
   values become missing and 2,463 observed Vs30 values change. The collection
   scripts were not updated between these uploads. Obtain the actual generating
   scripts and coordinate/Vs30 rules before claiming full reproducibility; the
   committed history does not establish exactly how those changes were made.

3. **Where should pipeline outputs be written and consumed?**
   Collection scripts write `data_processed`, while committed inputs live in
   `final_data`. Classification now reads `../final_data` from `pipeline_new`, but
   regression still calls `here("data_processed", ...)`. Agree on a canonical
   path and a documented execution directory. EDA already reads `final_data`.

## Model samples and comparisons

4. **Should both model families use one shared earthquake split?**
   Both currently have 106,877 eligible rows across 1,029 earthquakes before
   partitioning. Regression reconstructs an 823/206 event split (85,968/20,909
   rows); classification uses saved IDs, giving 760/191 represented events and
   79,084/23,348 rows. Another 4,445 eligible rows from 78 events are unassigned.
   There are 134 events in classification test/regression train, and another 134
   in classification train/regression test. This is not overlap within either
   method's own split, but their reported errors do not use the same test sample.
   Decide whether to share IDs or explicitly evaluate on common held-out events.
   Updating split IDs requires refitting and re-evaluating affected models.

5. **Can the committed classification figures be regenerated from the agreed run?**
   The saved Boosting Tier 1 confusion matrix totals 23,995 observations, whereas
   the current classification input/saved IDs yield 23,348. Refresh model figures
   from the approved data and split; do not relabel old figures as current output.
   The regression HTML has been updated to the raw-MMI model and its 823/206 split.

6. **Should Tier 3 use the same predictors as the other classification strategies?**
   `train_balanced` omits `report_count`, while the original `train` retains it.
   Unless this is intentional and documented, comparisons change both sampling
   and predictor availability and cannot isolate the effect of downsampling.

## Other exploratory entry point

7. **Should `data_exploration.Rmd` be retired or updated by its owner?**
   This separate legacy notebook reads the old path and refers to removed columns.
   It has not been rewritten as part of the DATA309_EDA update. The maintained EDA
   entry point and its generated products are documented in `DATA309_EDA.md`.

## Short message for the group

> I've updated DATA309_EDA to match the latest MMI >=3 preparation. Could we confirm
> the intended Mw threshold and share the scripts that generated the latest data,
> including the Vs30 changes? We also need to agree on the input path and a shared
> earthquake split before comparing models. The saved classification figures still
> use an older test sample, and Tier 3 currently drops report_count as a predictor.
