# Data-derived EDA text

Regenerate this text when inputs change; check the sample used by your models before copying it into the report.

The final dataset contains 473238 earthquake-cell observations. Requiring at least four reports per cell, complete selected fields, positive epicentral distance and Vs30, and MMI >=3 retains 106877 observations from 1029 earthquakes for regression.
In this sample, MMI 3 and 4 together account for 93.31% of observations, and 7154 observations have MMI 5 or higher. These are cell-level medians, not individual report ratings.
The final data contain 15353 missing Vs30 values. Classification uses MMI 3-8 and saved event partitions; 4445 otherwise eligible observations are outside those partitions. Identical eligible event-cell keys before splitting: TRUE.
The reconstructed regression split contains 85968 training observations across 823 earthquakes and 20909 test observations across 206 earthquakes. These partitions differ from the saved classification partitions.
The report-count retention analysis describes the sample-size trade-off; it does not establish that four reports is an optimal precision threshold.
Individual MMI frequencies are not included in the current inputs, so within-cell dispersion and comparisons of mean, median and mode cannot be estimated from these files.
