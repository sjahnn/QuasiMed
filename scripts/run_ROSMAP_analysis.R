################################################################################
# ROSMAP MedZIsc (Real Data) — pseudo-bulk feature construction + screening +
# final model + joint significance (Quasi / working-variance version)
#
# Author(s): Seungjun Ahn, Donald Porchia, Zhigang Li
#
#
# PURPOSE
#   This script applies the MedZIsc workflow to ROSMAP vasculature scRNA-seq data:
#     1) Build subject-level gene summaries from cell-level expression:
#        - M_g: subject-level mean expression among *nonzero* cells
#        - F_g: subject-level proportion of zeros (dropout rate), then logit(F_g)
#     2) Assemble an analysis-ready subject-level dataset:
#        - Outcome Y (here: PMI)
#        - Exposure X (here: AD diagnosis yes/no)
#        - Covariate(s) Z (here: sex)
#        - Mediators: {M_g, F_g} for filtered genes
#     3) Run MedZIsc:
#        - Outcome screening via LASSO: Y ~ [M, F]
#        - Mediator screening via quasi-GLM: M_g ~ X (+Z in JS stage),
#                                           F_g ~ X (+Z in JS stage)
#        - Combine-rule selection to cap features so final lm is estimable
#        - Final model: Y ~ X + Z + selected M + selected F
#        - Joint significance per mediator and BH FDR
#
# IMPORTANT NOTES FOR READERS
#   - Seurat `slot="data"` is typically normalized/log-transformed by upstream
#     preprocessing; this script uses the provided processed object as-is.
#   - cores/detectCores() are computed for SLURM/local environments, but this
#     particular script does not currently parallelize (foreach not used below).
#
# INPUTS
#   - Seurat object: Vasculature_cells.rds
#   - Subject metadata: individual_metadata_deidentified.tsv
#
# OUTPUT
#   - Saves `all_results` (list over parameter grid) as an .rds file in dir_results
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(MASS)     # NOTE: not used in the code below as written (kept for legacy)
  library(glmnet)   # used for cv.glmnet screening
  library(foreach)  # NOTE: not used in the code below as written (kept for legacy)
  library(dplyr)    # used for %>% arrange()/desc() in combine-rule ranking
})

# -------------------------------------------------------------------------
# Dataset provenance
# -------------------------------------------------------------------------
# Mathys et al. (Cell, 2023) ROSMAP single-cell atlas:
#   - cells annotated by subject and high-resolution cell type
#   - expression matrix stored in Seurat assays slot
# Supplementary: https://compbio.mit.edu/ad_aging_brain/#dataset-overlap-across-studies

################################################################################
# Section 0: Data inspection / preparation
################################################################################

# Cell-level metadata available in Seurat object:
#   - subject: unique (deidentified) ROSMAP participant ID for each cell
#   - cell_type_high_resolution: high-resolution cell type label
dat@meta.data[["subject"]]
dat@meta.data[["cell_type_high_resolution"]]

# Expression matrix:
#   NOTE: slot="data" usually corresponds to normalized/log-transformed values
#         in typical Seurat pipelines (depends on upstream processing).
counts <- GetAssayData(dat, slot = "data")
#dim(counts)  # genes x cells (e.g., 33538 x 17974)

################################################################################
# Main wrapper: run_singlecell_model()
#   Filters subjects/genes -> builds M/F -> runs MedZIsc -> returns results
################################################################################
run_singlecell_model <- function(min_Cells = 30,
                                 max_zeros  = 0.90,
                                 min_samples = 0.05) {

  # -------------------------------------------------------------------------
  # Filtering parameters (function arguments)
  #   min_Cells   : minimum number of cells required per subject
  #   max_zeros   : remove genes that are "too zero" for essentially all subjects
  #   min_samples : minimum fraction of subjects with expression > 0 (gene keep rule)
  # -------------------------------------------------------------------------

  # Re-load expression matrix inside function (uses global Seurat object `dat`)
  counts <- GetAssayData(dat, slot = "data")
  dim(counts)

  ################################################################################
  # Section 1: Build per-subject cell-by-gene matrices
  ################################################################################

  subjects      <- dat@meta.data[["subject"]]
  cell_barcodes <- colnames(counts)
  subject_ids   <- unique(subjects)

  # Create a named list of matrices:
  #   each element is a subject-specific matrix with rows=cells and columns=genes
  subject_expression_list <- lapply(subject_ids, function(subj) {
    idx <- which(subjects == subj)
    subj_cells <- cell_barcodes[idx]
    t(counts[, subj_cells, drop = FALSE])  # transpose: cell-by-gene
  })
  names(subject_expression_list) <- subject_ids

  # Remove subjects whose expression matrices sum to 0 (typically none after QC)
  subject_expression_list <- subject_expression_list[
    sapply(subject_expression_list, function(mat) sum(mat) > 0)
  ]

  # Diagnostic: list subjects removed (expected to be empty in this dataset)
  all_subjects     <- names(subject_expression_list)
  nonzero_subjects <- sapply(subject_expression_list, function(mat) sum(mat) > 0)
  removed_subjects <- all_subjects[!nonzero_subjects]
  print(removed_subjects)

  # Filter subjects by minimum cell count
  subject_expression_list_filtered <-
    subject_expression_list[sapply(subject_expression_list, nrow) >= min_Cells]

  # Number of retained subjects (post min_Cells filter)
  nSubjects <- length(names(subject_expression_list_filtered))


  ################################################################################
  # Section 2: Match subject-level metadata to expression list
  ################################################################################

  expr_ids <- names(subject_expression_list_filtered)
  meta_matched <- metadat[match(expr_ids, metadat$subject), ]
  rownames(meta_matched) <- expr_ids

  # Outcome Y: PMI
  outcome <- meta_matched$pmi

  # Exposure X: Pathologic diagnosis of AD (binary: yes=1, else=0)
  exposure_character <- meta_matched$Pathologic_diagnosis_of_AD
  exposure <- ifelse(exposure_character == "yes", 1, 0)

  # Covariate Z: sex (as provided; comment indicates 1=Male, 0=Female)
  covariate <- meta_matched$msex

  ################################################################################
  # Section 3: Derive subject-level mediator features (pseudo-bulk aggregation)
  #
  # For each subject i and gene g:
  #   M_{ig} = mean expression among nonzero cells (set 0 if all zeros)
  #   F_{ig} = proportion of zeros across cells, then clipped and logit-transformed
  #
  # NOTE:
  #   Because counts come from Seurat slot="data", M_{ig} is on whatever scale the
  #   upstream normalization produced (often log-normalized).
  ################################################################################

  # Subject-by-gene matrix of mean expression across all cells (not used downstream)
  M_subject_all <- t(sapply(subject_expression_list_filtered, function(x) Matrix::colMeans(x)))

  # Subject-by-gene matrix: mean across *nonzero* cells only
  M_nz_subject_all <- t(sapply(subject_expression_list_filtered, function(x) {
    col_sums   <- colSums(x * (x > 0))  # sum of nonzero values per gene
    col_counts <- colSums(x > 0)        # count of nonzero cells per gene
    ifelse(col_counts > 0, col_sums / col_counts, 0)
  }))

  # Subject-by-gene matrix of zero proportions
  F_subject_all <- t(sapply(subject_expression_list_filtered, function(x) Matrix::colMeans(x == 0)))

  # Clip F_g away from 0/1 to avoid infinite logit
  adjust_Fg <- function(Fg) {
    pmin(pmax(Fg, 0.001), 0.999)
  }
  F_subject_all <- apply(F_subject_all, 2, adjust_Fg)

  ################################################################################
  # Section 4: Gene filtering to reduce dimension and remove degenerate genes
  #
  # Filtering sequence:
  #   (A) Remove genes with total expression 0 across subjects
  #   (B) Remove genes that are "too zero" for all subjects (via max_zeros)
  #   (C) Keep genes expressed (>0) in at least min_samples fraction of subjects
  #
  #   Filtering is applied consistently to M and F to keep matched gene sets.
  ################################################################################

  # (A) Keep genes with nonzero total expression across subjects
  filtered_genes <- which(colSums(M_subject_all) > 0)
  M_subject_all     <- M_subject_all[, filtered_genes]
  M_nz_subject_all  <- M_nz_subject_all[, filtered_genes]
  F_subject_all     <- F_subject_all[, filtered_genes]

  # (B) Drop genes that have dropout proportion >= max_zeros for all subjects.
  #     Implementation: keep genes where at least one subject has F < max_zeros.
  filtered_genes <- which(apply(F_subject_all, 2, min) < max_zeros)
  M_subject_all    <- M_subject_all[, filtered_genes]
  M_nz_subject_all <- M_nz_subject_all[, filtered_genes]
  F_subject_all    <- F_subject_all[, filtered_genes]

  # (C) Keep genes expressed (>0) in > min_samples fraction of subjects
  # NOTE: This uses M_subject_all > 0 (all-cells mean), not M_nz_subject_all.
  filtered_genes <- which(colMeans(M_subject_all > 0) > min_samples)

  M_subject_all    <- M_subject_all[, filtered_genes]
  M_nz_subject_all <- M_nz_subject_all[, filtered_genes]
  F_subject_all    <- F_subject_all[, filtered_genes]

  # Final mediator matrices:
  M_subject <- M_nz_subject_all
  F_subject <- F_subject_all

  # Transform dropout proportion to logit scale: log(F/(1-F))
  F_subject <- log(F_subject / (1 - F_subject))

  ################################################################################
  # Section 5: Construct analysis-ready subject-level dataset
  ################################################################################

  # Prefix mediator columns to preserve gene identity and mediator type
  colnames(M_subject) <- paste0("M_", colnames(M_subject))
  colnames(F_subject) <- paste0("F_", colnames(F_subject))

  # Combine into a single dataframe
  #   Y: outcome (PMI)
  #   X: exposure (AD yes/no)
  #   Z: covariate (sex)
  rosmap_realdat <- data.frame(
    Y = outcome,
    X = exposure,
    Z = covariate,
    M_subject,
    F_subject
  )

  # Drop subjects with missing outcome
  na_rows <- which(is.na(rosmap_realdat[, "Y"]))
  rosmap_realdat <- rosmap_realdat[-na_rows, ]

  ################################################################################
  # Section 6: MedZIsc on ROSMAP data
  #
  # Screening + final model + joint significance testing, analogous to the
  # simulation evaluator but with:
  #   - covariate(s) Z included in final model and in X->M / X->F regressions
  #   - gene names are real gene symbols (after prefix stripping)
  ################################################################################

  set.seed(123)  # Ensures reproducible cv.glmnet folds
  n_genes <- ncol(rosmap_realdat[, grep("^(M_)", colnames(rosmap_realdat))])

  MedZIsc <- function(data.name, n_genes, covariate.names) {

    # Design matrix for LASSO screening: all M_ and F_ columns
    MF_mat <- as.matrix(data.name[, grep("^(M_|F_)", colnames(data.name))])
    Y_vec  <- data.name$Y

    # -----------------------------------------------------------------------
    # STEP 1: Outcome screening via LASSO
    #   Fit: Y ~ [M,F] (no Z here; Z comes in the final model)
    # -----------------------------------------------------------------------
    outcome.lasso <- cv.glmnet(MF_mat, Y_vec, alpha = 1, type.measure = "mse")
    coef_lasso <- coef(outcome.lasso, s = "lambda.min")

    selected_genes <- as.matrix(coef_lasso)
    selected_genes <- selected_genes[selected_genes[, 1] != 0, , drop = FALSE]
    selected_gene_names <- setdiff(rownames(selected_genes), "(Intercept)")

    # -----------------------------------------------------------------------
    # STEP 2: Mediator screening via per-gene quasi GLMs
    #   M-side: M_g ~ X
    #   F-side: F_g ~ X
    #   NOTE: working variance is constant, identity link (treating as continuous).
    # -----------------------------------------------------------------------

    # ---- M side screening ----
    lm_models_M <- list()
    lm_coeff_M  <- numeric(n_genes)
    pvals_M     <- numeric(n_genes)

    selected_M_gene_names <- character(0)
    selected_M_gene_pvals <- numeric(0)

    for (g in 1:n_genes) {
      gene_col <- colnames(data.name[, grep("^(M_)", colnames(data.name))])[g]
      formula_str <- as.formula(paste(gene_col, "~ X"))

      fit <- suppressWarnings({
        tryCatch(
          glm(formula_str, family = quasi(link = "identity", variance = "constant"), data = data.name),
          error = function(e) NULL
        )
      })

      if (!is.null(fit)) {
        lm_models_M[[g]] <- fit
        lm_coeff_M[g] <- coef(fit)["X"]
        pvals_M[g] <- summary(fit)$coefficients["X", "Pr(>|t|)"]

        if (!is.na(pvals_M[g]) && pvals_M[g] < 0.05) {
          selected_M_gene_names <- c(selected_M_gene_names, gene_col)
          selected_M_gene_pvals <- c(selected_M_gene_pvals, pvals_M[g])
        }
      } else {
        lm_models_M[[g]] <- NULL
        lm_coeff_M[g] <- NA
        pvals_M[g] <- NA
      }
    }

    M_gene_pval_table <- data.frame(pval = selected_M_gene_pvals)
    rownames(M_gene_pval_table) <- selected_M_gene_names

    # ---- F side screening ----
    beta_models_F <- list()
    beta_coeff_F  <- numeric(n_genes)
    pvals_F       <- numeric(n_genes)

    selected_F_gene_names <- character(0)
    selected_F_gene_pvals <- numeric(0)

    for (g in 1:n_genes) {
      gene_col <- colnames(data.name[, grep("^(F_)", colnames(data.name))])[g]
      formula_str <- as.formula(paste(gene_col, "~ X"))

      # capture.output/suppressMessages is used to silence noisy model output
      invisible(
        suppressWarnings(
          suppressMessages(
            capture.output({
              fit <- tryCatch(
                glm(formula_str, family = quasi(link = "identity", variance = "constant"), data = data.name),
                error = function(e) NULL
              )
            }, type = "message")
          )
        )
      )

      if (!is.null(fit)) {
        beta_models_F[[g]] <- fit
        beta_coeff_F[g] <- coef(fit)["X"]
        pvals_F[g] <- summary(fit)$coefficients["X", "Pr(>|t|)"]

        if (!is.na(pvals_F[g]) && pvals_F[g] < 0.05) {
          selected_F_gene_names <- c(selected_F_gene_names, gene_col)
          selected_F_gene_pvals <- c(selected_F_gene_pvals, pvals_F[g])
        }
      } else {
        beta_models_F[[g]] <- NULL
        beta_coeff_F[g] <- NA
        pvals_F[g] <- NA
      }
    }

    F_gene_pval_table <- data.frame(pval = selected_F_gene_pvals)
    rownames(F_gene_pval_table) <- selected_F_gene_names

    # -----------------------------------------------------------------------
    # STEP 3: Combine rule (cap feature count so final lm is estimable)
    #
    # Convert selected variable names to gene identifiers (strip "M_" / "F_")
    # -----------------------------------------------------------------------
    Y_indices <- unique(sub("^[^_]*_", "", selected_gene_names))
    M_indices <- sub("^[^_]*_", "", selected_M_gene_names)
    F_indices <- sub("^[^_]*_", "", selected_F_gene_names)

    # Number of mediators to keep per side:
    #   approx n/4 - (#covariates) - 1 (exposure)  [legacy rule]
    num <- round(nrow(data.name) / 4 - length(covariate.names) - 1)

    # -------------------- M-side combine/rank --------------------
    # Priority:
    #   1) Y ∩ M (rank by LASSO coefficient)
    #   2) Y \ M (rank by LASSO coefficient)
    #   3) M \ Y (rank by smallest p-value)
    Y_intersect_M <- intersect(Y_indices, M_indices)

    # NOTE: `formatching` is used as a temporary numeric vector for ordering.
    formatching <- selected_genes[
      which(sub("M_", "", rownames(selected_genes)) %in% Y_intersect_M),
    ]
    YM_model_coeff <- data.frame(
      formatching,
      row.names = rownames(selected_genes)[
        which(sub("M_", "", rownames(selected_genes)) %in% Y_intersect_M)
      ]
    )
    YM_model_coeff <- data.frame(YM_model_coeff %>% arrange(desc(formatching)))
    rownames(YM_model_coeff) <- sub("M_", "", rownames(YM_model_coeff))

    Y_butnot_M <- setdiff(Y_indices, M_indices)
    Y_butnot_M_coeff <- data.frame(
      formatching = selected_genes[
        which(sub("M_", "", rownames(selected_genes)) %in% Y_butnot_M),
      ],
      row.names = rownames(selected_genes)[
        which(sub("M_", "", rownames(selected_genes)) %in% Y_butnot_M)
      ]
    )
    Y_butnot_M_coeff <- data.frame(Y_butnot_M_coeff %>% arrange(desc(formatching)))
    rownames(Y_butnot_M_coeff) <- sub("M_", "", rownames(Y_butnot_M_coeff))

    M_butnot_Y <- setdiff(M_indices, Y_indices)
    formatching <- M_gene_pval_table[
      which(sub("M_", "", rownames(M_gene_pval_table)) %in% M_butnot_Y),
    ]
    M_butnot_Y_pval <- data.frame(formatching)
    rownames(M_butnot_Y_pval) <- M_butnot_Y
    M_butnot_Y_pval <- data.frame(M_butnot_Y_pval %>% arrange(formatching))

    combined_tab_M <- rbind(YM_model_coeff, Y_butnot_M_coeff, M_butnot_Y_pval)
    M_names <- rownames(combined_tab_M)[1:min(num, nrow(combined_tab_M))]

    # -------------------- F-side combine/rank --------------------
    Y_intersect_F <- intersect(Y_indices, F_indices)

    YF_model_coeff <- data.frame(
      formatching = selected_genes[
        which(sub("F_", "", rownames(selected_genes)) %in% Y_intersect_F),
      ],
      row.names = rownames(selected_genes)[
        which(sub("F_", "", rownames(selected_genes)) %in% Y_intersect_F)
      ]
    )
    YF_model_coeff <- data.frame(YF_model_coeff %>% arrange(desc(formatching)))
    rownames(YF_model_coeff) <- sub("F_", "", rownames(YF_model_coeff))

    Y_butnot_F <- setdiff(Y_indices, F_indices)
    Y_butnot_F_coeff <- data.frame(
      formatching = selected_genes[
        which(sub("F_", "", rownames(selected_genes)) %in% Y_butnot_F),
      ],
      row.names = rownames(selected_genes)[
        which(sub("F_", "", rownames(selected_genes)) %in% Y_butnot_F)
      ]
    )
    Y_butnot_F_coeff <- data.frame(Y_butnot_F_coeff %>% arrange(desc(formatching)))
    rownames(Y_butnot_F_coeff) <- sub("F_", "", rownames(Y_butnot_F_coeff))

    F_butnot_Y <- setdiff(F_indices, Y_indices)
    F_butnot_Y_pval <- data.frame(
      formatching = F_gene_pval_table[
        which(sub("F_", "", rownames(F_gene_pval_table)) %in% F_butnot_Y),
      ]
    )
    rownames(F_butnot_Y_pval) <- F_butnot_Y
    F_butnot_Y_pval <- data.frame(F_butnot_Y_pval %>% arrange(formatching))

    combined_tab_F <- rbind(YF_model_coeff, Y_butnot_F_coeff, F_butnot_Y_pval)
    F_names <- rownames(combined_tab_F)[1:min(num, nrow(combined_tab_F))]

    # -----------------------------------------------------------------------
    # STEP 4: Final model fit
    #   Y ~ X + Z + selected M + selected F
    # -----------------------------------------------------------------------
    Z_names <- covariate.names
    M_names <- if (length(M_names) > 0) paste0("M_", M_names) else character(0)
    F_names <- if (length(F_names) > 0) paste0("F_", F_names) else character(0)
    predictors <- c("X", Z_names, M_names, F_names)

    if (length(M_names) > 0 | length(F_names) > 0) {

      final_formula <- as.formula(paste("Y ~", paste(predictors, collapse = " + ")))
      fit_final_model <- lm(final_formula, data = data.name)

      coef_outcome <- coef(fit_final_model)
      se_outcome   <- coef(summary(fit_final_model))[, "Std. Error"]

      # NOTE: Uses normal approximation p-values rather than t-distribution p-values
      pvals_final <- 2 * (1 - pnorm(abs(coef_outcome / se_outcome)))

      # ---------------------------------------------------------------------
      # STEP 5: Joint Significance Testing (JS)
      #
      # For each selected mediator:
      #   M-side: p_g^M = max( p(X->M_g | Z), p(M_g->Y | X,Z,selected set) )
      #   F-side: p_g^F = max( p(X->F_g | Z), p(F_g->Y | X,Z,selected set) )
      #
      # BH adjustment is applied over the selected mediator sets (as implemented).
      #
      # NOTE: The "NaN" handling below compares to the string "NaN".
      #       In base R, NaN is typically numeric; a safer check is is.nan().
      #       Kept unchanged here to preserve existing behavior.
      # ---------------------------------------------------------------------
      joint_pvals_M <- joint_pvals_F <- gamma_coef_vec <- alpha_coef_vec <- NULL
      gamma_se_vec  <- alpha_se_vec  <- NULL

      for (m_var in M_names) {
        formula_Mg <- as.formula(paste0(m_var, " ~ X + ", paste(Z_names, collapse = " + ")))
        fit_Mg <- glm(formula_Mg, family = quasi(link = "identity", variance = "constant"), data = data.name)

        gamma_coef <- coef(fit_Mg)["X"]
        p_gamma_Xg <- ifelse(
          summary(fit_Mg)$coefficients["X", "Pr(>|t|)"] == "NaN",
          1,
          summary(fit_Mg)$coefficients["X", "Pr(>|t|)"]
        )

        gamma_coef_vec <- c(gamma_coef_vec, gamma_coef)
        gamma_se <- coef(summary(fit_Mg))["X", "Std. Error"]
        gamma_se_vec <- c(gamma_se_vec, gamma_se)

        p_beta_Mg <- pvals_final[m_var]
        joint_pvals_M <- c(joint_pvals_M, max(p_gamma_Xg, p_beta_Mg, na.rm = TRUE))
      }

      for (f_var in F_names) {
        formula_Fg <- as.formula(paste0(f_var, " ~ X + ", paste(Z_names, collapse = " + ")))
        fit_Fg <- glm(formula_Fg, family = quasi(link = "identity", variance = "constant"), data = data.name)

        p_alpha_Xg <- ifelse(
          summary(fit_Fg)$coefficients["X", "Pr(>|t|)"] == "NaN",
          1,
          summary(fit_Fg)$coefficients["X", "Pr(>|t|)"]
        )

        alpha_coef <- coef(fit_Fg)["X"]
        alpha_coef_vec <- c(alpha_coef_vec, alpha_coef)
        alpha_se <- coef(summary(fit_Fg))["X", "Std. Error"]
        alpha_se_vec <- c(alpha_se_vec, alpha_se)

        p_beta_Fg <- pvals_final[f_var]
        joint_pvals_F <- c(joint_pvals_F, max(p_alpha_Xg, p_beta_Fg, na.rm = TRUE))
      }

      if (!is.null(joint_pvals_M)) names(joint_pvals_M) <- M_names
      if (!is.null(joint_pvals_F)) names(joint_pvals_F) <- F_names

      # BH FDR control over joint p-values
      fdr_M <- p.adjust(joint_pvals_M, method = "BH")
      fdr_F <- p.adjust(joint_pvals_F, method = "BH")
      signif_M <- which(fdr_M < 0.05)
      signif_F <- which(fdr_F < 0.05)

      # Return structure varies depending on whether M-side/F-side mediators exist
      if (is.null(joint_pvals_M) & !is.null(joint_pvals_F)) {
        return(list(
          nSubjects     = nSubjects,
          n_genes       = n_genes,
          coef_outcome  = coef_outcome,
          se_outcome    = se_outcome,
          alpha_coef_vec = alpha_coef_vec,
          alpha_se_vec   = alpha_se_vec,
          signif_F      = signif_F,
          adj_pvals_F   = fdr_F,
          joint_pvals_F = joint_pvals_F
        ))
      }

      if (!is.null(joint_pvals_M) & is.null(joint_pvals_F)) {
        return(list(
          nSubjects     = nSubjects,
          n_genes       = n_genes,
          coef_outcome  = coef_outcome,
          se_outcome    = se_outcome,
          gamma_coef_vec = gamma_coef_vec,
          gamma_se_vec   = gamma_se_vec,
          signif_M      = signif_M,
          adj_pvals_M   = fdr_M,
          joint_pvals_M = joint_pvals_M
        ))
      }

      if (!is.null(joint_pvals_M) & !is.null(joint_pvals_F)) {
        return(list(
          nSubjects      = nSubjects,
          n_genes        = n_genes,
          coef_outcome   = coef_outcome,
          se_outcome     = se_outcome,
          gamma_coef_vec = gamma_coef_vec,
          gamma_se_vec   = gamma_se_vec,
          alpha_coef_vec = alpha_coef_vec,
          alpha_se_vec   = alpha_se_vec,
          signif_M       = signif_M,
          signif_F       = signif_F,
          adj_pvals_M    = fdr_M,
          adj_pvals_F    = fdr_F,
          joint_pvals_M  = joint_pvals_M,
          joint_pvals_F  = joint_pvals_F
        ))
      }

    } else {
      message("No mediator detected.")
      return(list(
        nSubjects      = nSubjects,
        n_genes        = n_genes,
        coef_outcome   = NULL,
        se_outcome     = NULL,
        gamma_coef_vec = NULL,
        gamma_se_vec   = NULL,
        alpha_coef_vec = NULL,
        alpha_se_vec   = NULL,
        signif_M       = integer(0),
        signif_F       = integer(0),
        adj_pvals_M    = NULL,
        adj_pvals_F    = NULL,
        joint_pvals_M  = NULL,
        joint_pvals_F  = NULL
      ))
    }
  }

  # Run MedZIsc once for this filtered/constructed dataset
  ROSMAPresults <- MedZIsc(data.name = rosmap_realdat, n_genes = n_genes, covariate.names = "Z")

  # Tag filtering parameters into the result for traceability
  ROSMAPresults$min_Cells   <- min_Cells
  ROSMAPresults$max_zeros   <- max_zeros
  ROSMAPresults$min_samples <- min_samples

  return(ROSMAPresults)
}

################################################################################
# Parameter grid runner
#   Runs run_singlecell_model() over a grid of filtering parameters
#   (here, the grid is effectively a single setting)
################################################################################
param_grid <- expand.grid(
  min_Cells   = seq(30, 30, by = 10),
  max_zeros   = seq(0.9, 0.9, by = 0.05),
  min_samples = seq(0.05, 0.05, by = 0.05),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

all_results <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  pars <- param_grid[i, ]

  all_results[[i]] <- run_singlecell_model(
    min_Cells   = pars$min_Cells,
    max_zeros   = pars$max_zeros,
    min_samples = pars$min_samples
  )

  # Convenience: sort adjusted p-values and joint p-values for easy inspection
  all_results[[i]]$adj_pvals_M   <- sort(all_results[[i]]$adj_pvals_M)
  all_results[[i]]$adj_pvals_F   <- sort(all_results[[i]]$adj_pvals_F)
  all_results[[i]]$joint_pvals_M <- sort(all_results[[i]]$joint_pvals_M)
  all_results[[i]]$joint_pvals_F <- sort(all_results[[i]]$joint_pvals_F)
}

################################################################################
# Save results
################################################################################
cur_time <- format(Sys.time(), "%Y_%m_%d_%H%M")
saveRDS(
  all_results,
  file = file.path(dir_results, paste0("quasimed_Rosmap_Results", "_", cur_time, ".rds"))
)

