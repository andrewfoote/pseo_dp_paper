# Description of Programs

This folder includes all the programs used to run the simulations described in the paper. They are numbered sequentially on purpose.

## 00.truepercentiles.sas

This program takes the input data, and calculates the true percentiles at the cell level. We use these as the baseline for what the "truth" is, when we construct our error measures.

## 01.iteration_percentiles.sas

This was used earlier in the paper- no longer in use.

## 02.iteration_percentiles.sas

This program iterates over three loops for two different types of histograms, log normal and evenly spaced bins. We iterate over all the combinations of the following: numbins x epsilon x iteration (20 total iterations for each numbin x epsilon combination).

We then save those results for later, by taking the mean of relative and L1 errors, and relative accuracy, over the following values: iteration x numbin x epsilon x cellbin, where cellbin distinguishes different cell sizes.

## 02.1.ssiteration_percentiles.sas

This program calculates smooth sensitivity for each cell, then takes 20 draws from the the distribution 1/|1+x|^4, to calculate the protected percentile value. At end, we calculate different error and accuracy measures from comparison purposes.

## 03.percentile_errors_export.sas

This program exports the errors measures into STATA DTA files. It is not necessary, but I prefer the graphing tools in STATA.
