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
-  `run_singlecell_model()`: This a wrapper function for applying the MedZIsc workflow to single-cell data after subject-level aggregation. 
It takes cell-level expression data, filters subjects and genes using user-specified thresholds, constructs two gene-level mediator features for each subject, 
and prepares the final analysis dataset. For each gene, the function derives an **M** feature (mean expression among nonzero cells) and an **F** feature (proportion of zero-expression cells),
where **F** is transformed using the logit function. The resulting subject-level dataset is then passed to `run_quasimed()`, and the fitted results are returned together with the filtering settings used.

-  `run_quasimed()`: This is the main function implementing quasi-regression mediation analysis with paired mediators derived from zero-inflated single-cell data. It first performs outcome screening using LASSO,
  then screens the mediator and zero-inflation components separately through gene-wise quasi-regression models. The final model includes the exposure, covariates, and the selected mediator features.
Joint significance tests are then performed for the **M** and **F** components of each selected gene, and the function returns coefficient estimates, standard errors, and Benjamini-Hochberg adjusted p-values.

## Important Notes
-  For a step-by-step example of the full analysis workflow, see `scripts/run_rosmap_analysis.R`. This script illustrates how the single-cell data are processed using **Seurat** R package, 
how subject-level mediator features are constructed, and how additional filtering steps are applied before running the QuasiMed (or `run_quasimed`). The final ROSMAP data analysis can be performed either 
directly with `run_quasimed()` or through the wrapper function `run_singlecell_model()`.

- We used pre-processed single-cell gene count data for vascular and epithelial cells (Vasculature cells.rds) together with de-identified clinical metadata (individual metadata deidentified.tsv) from the ROSMAP study.
- Please download these data from the AD/Aging Brain Atlas (https://compbio.mit.edu/ad_aging_brain/). We do not own the data, so the data is not saved in this GitHub repository.
