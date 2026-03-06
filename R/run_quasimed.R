#' @title run_quasimed
#'
#' @description A main function for conducting causal mediation analysis
#' under the quasi-regression mediation framework for zero-inflated
#' single-cell data.
#' @param data.name A data.frame or matrix with N x (2G + K), where N
#'   is the number of samples, G is the number of genes (each gene contributes
#'   two features: one for the zero component and one for the non-zero
#'   component), and K is the number of covariates.
#' @param n_genes An integer value. The number of genes (G) represented
#'   in the data.
#' @param covariate.names A character vector to specify the column name of covariates.
#' @return A list containing: (1) number of genes; (2) estimated coefficients from
#' the outcome and two mediation models (M and F models in the methodology paper);
#' (3) standard errors corresponding to the estimated coefficients;
#' (4) logical vector indicating whether each gene's mediator
#' component (M model) is statistically significant;
#' (5) logical vector indicating whether each gene's zero-inflation
#' component (F model) is statistically significant;
#' (6) adjusted p-values for M and F models (joint significance test).
#' @examples
#' \donttest{
#'   set.seed(123)  # Ensures reproducible cv.glmnet folds
#'   n_genes <- ncol(rosmap_realdat[, grep("^(M_)", colnames(rosmap_realdat))])
#'   run_quasimed(data.name = simulated_data, n_genes = n_genes, covariate.names = c("Z1", "Z2", "Z3"))
#' }
#' @import glmnet
#' @importFrom stats as.formula coef lm p.adjust pnorm
#' @importFrom utils capture.output
#' @export

run_quasimed <- function(data.name, n_genes, covariate.names) {

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
