# DATA309 EDA outputs

Generated: 2026-10-17 08:03:13 UTC

Source: current final_data Parquet files. Observation: one earthquake-cell pair; response: median reported MMI.

## Current sample

Final data: 473238 earthquake-cells across 1584 earthquakes, from 2016-09-01 16:37:55 to 2026-10-15 05:25:34 (UTC).
Regression preparation retains 106878 complete earthquake-cells across 1029 earthquakes with at least four reports per cell.
Classification preparation additionally requires MMI levels 3-8: 106877 eligible cells; 4445 are not covered by the saved IDs.
The two samples and partitions must not be described as identical. See 01_split_comparison.csv and 13_pipeline_audit.csv.

## Reading the results

00: sample flow and catalogue coverage. 01: saved-ID coverage and reconstructed split comparison.
02-06: MMI, missingness, threshold retention, event coverage and distributions.
07-10: descriptive distance/magnitude/site associations and within-event examples.
11: optional local coverage map. 12: unavailable analyses. 13: model integration issues outside EDA scope.
All percentages in 02 use their own stage as denominator; 04 uses all final earthquake-cells.
Vs30 and response summaries count earthquake-cells, so repeated locations/events receive repeated weight.
Cell IDs in the current inputs include the earthquake identifier; distinct_cell_ids is not a count of unique physical locations. reporting_coordinates counts distinct coordinate pairs.
Spearman correlations require at least ten observations and variation in both variables. NA means not estimable, not zero.
Full-sample EDA is descriptive; it is not held-out performance validation or evidence for tuning against test observations.

## Inputs and reproducibility

input_manifest.csv identifies input and source-code hashes. output_manifest.csv identifies this run's products.
No models are fitted, no data are downloaded, and no input datasets or saved IDs are overwritten.
Both model notebooks still refer to an old data_processed directory; their execution is outside this EDA update.
See ../DATA309_EDA.md for commands, dependencies, migration and limitations.
