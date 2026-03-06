# QuasiMed

## Overview:
A mediation framework specialized for single-cell data that comprises the following: 
(1) mediator candidates screening through penalization and marginal modeling, (2) estimation of indirect effects through the average expression and proportion of expressing cells, 
and (3) hypothesis testing with multiciplity control. Full methodological details are available in our recent preprint.

## How to Install this R Package:
```r
## From GitHub: 
devtools::install_github("sjahnn/QuasiMed")
```

## Key Functions
-  `run_singlecell_model()`: This a wrapper function for applying the QuasiMed workflow to single-cell data after aggregating to the subject level. It filters subjects and genes based on user-defined thresholds, then builds the subject-level inputs used for mediation analysis. For each gene, the function computes two mediator variables per subject: **M**, the mean expression across cells with nonzero expression, and **F**, the fraction of cells with zero expression. **F** is logit-transformed before analysis. The aggregated dataset is then passed to `run_quasimed()`, and the function returns both the fitted results and the filtering parameters that were applied.

-  `run_quasimed()`: This is the main function implementing quasi-regression mediation analysis with paired mediators derived from zero-inflated single-cell data. It first performs preliminary mediator screening through LASSO penalization (outcome model) and separately through gene-wise quasi-regression models (marginal modeling). The final model includes the exposure, covariates, and the selected mediator features. Joint significance tests are then performed for the **M** and **F** components of each selected gene, and the function returns coefficient estimates, standard errors, and Benjamini-Hochberg adjusted p-values.

## Important Notes
-  For a step-by-step example of the full analysis workflow, see `scripts/run_rosmap_analysis.R`. This script details how the single-cell data are processed using **Seurat** R package, 
how subject-level mediator features are constructed, and how additional filtering steps are applied before running the QuasiMed (or `run_quasimed`). The final ROSMAP data analysis can be performed either 
directly with `run_quasimed()` or through the wrapper function `run_singlecell_model()`.

- We used pre-processed single-cell gene count data for vascular and epithelial cells (Vasculature cells.rds) together with de-identified clinical metadata (individual metadata deidentified.tsv) from the ROSMAP study.
- Please download these data from the AD/Aging Brain Atlas (https://compbio.mit.edu/ad_aging_brain/). We do not own the data, so the data is not saved in this GitHub repository.
