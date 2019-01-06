# Resampling-Tests-for-Lasso-Supplementary Material

This repository contains R code and tutorials for fitting the LassoPL model from Arbet et al "Resampling-based tests for Lasso in
genome-wide association studies" BMC Genetics.  Accepted June 2017.

For small to moderate datasets (e.g. less than 100,000 predictors), refer to glmnet_lassoPL_functions.R and glmnet_lassoPL_tutorial.R.  These files use the glmnet R package and require loading the entire dataset into RAM.  Parallel computing is not used in these files.

For larger datasets(e.g. more than 100,000 predictors), refer to biglassoPL_functions.R and biglassoPL_tutorial.R.  These files use the bigmemory and biglasso packages to avoid loading the entire dataset into RAM.  Parallel computing is supported.

Any questions can be emailed to Jaron Arbet at: jaron.arbet@ucdenver.edu
