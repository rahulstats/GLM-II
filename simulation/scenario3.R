library(glmnet)
library(Matrix)
library(parallel)

# Clean workspace
graphics.off()  
rm(list = ls()) 

#-----------------------------#
# 1. Ultra-Fast Sparse Design Matrices
#-----------------------------#
L_mat_sparse <- function(k, m) {
  kronecker(Diagonal(choose(k, 2)), matrix(1, nrow = m, ncol = 1))
}

B_mat_sparse <- function(K) {
  pairs <- t(combn(K, 2))
  K2 <- nrow(pairs)
  sparseMatrix(
    i = rep(1:K2, 2),
    j = c(pairs[, 1], pairs[, 2]),
    x = c(rep(1, K2), rep(-1, K2)),
    dims = c(K2, K)
  )
}

C_mat_sparse <- function(K) {
  pairs <- t(combn(K, 2))
  triads <- t(combn(K, 3))
  K2 <- nrow(pairs)
  K3 <- nrow(triads)
  
  pair_mat <- matrix(0, K, K)
  pair_mat[lower.tri(pair_mat)] <- 1:K2
  pair_mat <- pair_mat + t(pair_mat)
  
  e1 <- pair_mat[cbind(triads[,1], triads[,2])]
  e2 <- pair_mat[cbind(triads[,1], triads[,3])]
  e3 <- pair_mat[cbind(triads[,2], triads[,3])]
  
  sparseMatrix(
    i = c(e1, e2, e3),
    j = rep(1:K3, 3),
    x = c(rep(1, K3), rep(-1, K3), rep(1, K3)),
    dims = c(K2, K3)
  )
}

mu_par <- function(K) { combn(K, 2, FUN = function(x) x[2] - x[1]) }

#-----------------------------#
# 2. Parameter Initialization
#-----------------------------#
K = 20   
m = 30   
rep_sim = 1000 

K2 = choose(K, 2)
K3 = choose(K, 3)

cat("Initializing matrices for K =", K, "... (Edges:", K2, "| Triads:", K3, ")\n")

BK = B_mat_sparse(K)
CK = C_mat_sparse(K)
L_Mat = L_mat_sparse(K, m)
LC = L_Mat %*% CK 
colnames(CK) = paste("v", 1:ncol(CK), sep="")

#-----------------------------#
# EXACT LEXICOGRAPHICAL INDEXING
#-----------------------------#
# True Model: 1*(1,2,3) - 1*(1,2,4) + 1*(1,2,5)
a = 1; b = -1

# Dynamically find the exact column indices for (1,2,3) and (1,2,4)
all_triads <- combn(K, 3)
idx1 <- 1; idx2 <- 2; idx3 <- 3

true_triads <- c(idx1, idx2, idx3) 
n_true <- length(true_triads) 

cat("True Triad Indices mapping to (1,2,3), (1,2,4) and (1,2,5):", true_triads, "\n")

nu = as.vector(mu_par(K)/10 + a*CK[,true_triads[1]] + b*CK[,true_triads[2]]) + a*CK[,true_triads[3]]
Y_signal = as.vector(L_Mat %*% nu)

#-----------------------------#
# 3. Fast Precomputations
#-----------------------------#
BK_mat <- as.matrix(BK)
BK_full_rank <- BK_mat[, -K] 
tBK_BK_inv <- chol2inv(chol(crossprod(BK_full_rank)))
P_BK <- BK_full_rank %*% tBK_BK_inv %*% t(BK_full_rank)

var_cyc_nu <- rep(K2, K2) 
Pair_in_Triad_mat <- abs(CK)
c_xmat <- 1 / K2 
CK_t <- as.matrix(CK[ ,true_triads])

#-----------------------------#
# 4. Multi-Core Setup
#-----------------------------#
num_cores <- min(8, max(1, detectCores() - 1)) 
cat("Ready. Using", num_cores, "cores. Starting Simulation...\n")

RNGkind("L'Ecuyer-CMRG")
set.seed(2345)

#-----------------------------#
# 5. Parallel Simulation Loop
#-----------------------------#
results_list <- mclapply(1:rep_sim, function(ai) {
  
  tryCatch({
    
    out_n_triad_lasso <- NA; out_true_model_select_lasso <- NA
    out_n_triad_hybrid <- NA; out_true_model_select_hybrid <- NA 
    out_n_triad_ftbs <- NA; out_true_model_select_ftbs <- NA 
    
    bY <- as.vector(Y_signal + rnorm(K2 * m, 0, 1))
    nu_hat <- as.vector(crossprod(L_Mat, bY) / m)
    
    hat_nu_linear <- as.vector(P_BK %*% nu_hat)
    hat_nu_cyc <- nu_hat - hat_nu_linear 
    
    #-------------------------------------------------#
    # EXACT LINEAR SPAN CHECKER
    #-------------------------------------------------#
    check_span <- function(cols) {
      if(length(cols) == 0) return(0)
      X_eval <- cbind(BK_full_rank, as.matrix(CK[, cols, drop=FALSE]))
      fit_eval <- lm.fit(X_eval, nu) 
      coefs_eval <- ifelse(is.na(fit_eval$coefficients), 0, fit_eval$coefficients)
      res_ss <- sum((nu - as.vector(X_eval %*% coefs_eval))^2)
      return(ifelse(res_ss < 1e-8, n_true, 0)) 
    }
    
    #=================================================#
    # METHOD 1: STANDARD LASSO
    #=================================================#
    bY_adj <- as.vector(bY - L_Mat %*% hat_nu_linear)
    cv_model <- cv.glmnet(LC, bY_adj, alpha = 1, intercept = F, nfolds = 3) 
    best_model <- glmnet(LC, bY_adj, alpha = 1, lambda = cv_model$lambda.min, intercept = F)
    
    lasso_coefs <- coef(best_model)[-1]
    lasso_select <- which(lasso_coefs != 0)
    lasso_select_values <- abs(lasso_coefs[lasso_select])
    
    out_n_triad_lasso <- length(lasso_select)
    out_true_model_select_lasso <- check_span(lasso_select)
    
    #=================================================#
    # METHOD 2: HYBRID (Post-LASSO Pruning)
    #=================================================#
    if (length(lasso_select) > 0) {
      mat_x_hybrid <- as.matrix(CK[, lasso_select, drop=FALSE])
      qrx_h <- qr(mat_x_hybrid)
      pivot_idx_h <- qrx_h$pivot[seq_len(qrx_h$rank)]
      
      final_MatX_hybrid <- mat_x_hybrid[, pivot_idx_h, drop=FALSE]
      current_cols_hybrid <- lasso_select[pivot_idx_h]
      current_lasso_vals <- lasso_select_values[pivot_idx_h]
      
      keep_idx_hybrid <- 1:length(current_cols_hybrid)
      
      while(length(keep_idx_hybrid) > 0) {
        BC_X <- cbind(BK_full_rank, final_MatX_hybrid[, keep_idx_hybrid, drop=FALSE])
        p_X <- ncol(BC_X)
        df_resid <- K2 - p_X
        
        if (df_resid <= 0) {
          worst_idx <- which.min(current_lasso_vals[keep_idx_hybrid])
          keep_idx_hybrid <- keep_idx_hybrid[-worst_idx]
          next
        }
        
        fit_X <- lm.fit(BC_X, nu_hat)
        
        # Rank Deficiency Check
        triad_coefs <- fit_X$coefficients[K:p_X]
        if (any(is.na(triad_coefs))) {
          worst_idx <- which(is.na(triad_coefs))[1]
          keep_idx_hybrid <- keep_idx_hybrid[-worst_idx]
          next
        }
        
        sig2 <- sum(fit_X$residuals^2) / df_resid
        XtX_inv <- chol2inv(chol(crossprod(BC_X)))
        se <- sqrt(sig2 * diag(XtX_inv))
        
        triad_idx <- K:p_X
        t_stats <- fit_X$coefficients[triad_idx] / se[triad_idx]
        pvals <- 2 * pt(-abs(t_stats), df = df_resid)
        
        if (max(pvals) > 0.01) {
          keep_idx_hybrid <- keep_idx_hybrid[-which.max(pvals)]
        } else {
          break
        }
      }
      
      if (length(keep_idx_hybrid) > 0) {
        final_cols_h <- current_cols_hybrid[keep_idx_hybrid]
        out_n_triad_hybrid <- length(final_cols_h)
        out_true_model_select_hybrid <- check_span(final_cols_h)
      } else {
        out_n_triad_hybrid <- 0
        out_true_model_select_hybrid <- 0
      }
    } else {
      out_n_triad_hybrid <- 0
      out_true_model_select_hybrid <- 0
    }
    
    #=================================================#
    # METHOD 3: PURE FTBS (FDR + Pruning)
    #=================================================#
    abs_t_stats_cyc <- abs(sqrt(m * K2) * hat_nu_cyc / sqrt(var_cyc_nu))
    nu_cyc_pval <- pnorm(abs_t_stats_cyc, mean=0, sd=1, lower.tail = F)
    fdr_pvals <- p.adjust(nu_cyc_pval, method = "BH")
    sig_flags <- as.numeric(fdr_pvals < 0.05) 
    
    tick_tab_colsum <- as.vector(crossprod(Pair_in_Triad_mat, sig_flags))
    tick3_col <- which(tick_tab_colsum == 3)
    tick2_col <- which(tick_tab_colsum == 2)
    tick1_col <- which(tick_tab_colsum == 1)
    
    triad_scores <- as.vector(crossprod(Pair_in_Triad_mat, abs_t_stats_cyc))
    
    test_model_fit <- function(cols_to_use) {
      if (length(cols_to_use) == 0 || length(cols_to_use) >= K2) {
        return(list(pval = 0, MatX = NULL, cols = NULL))
      }
      
      mat_x <- as.matrix(CK[, cols_to_use, drop=FALSE])
      if(length(cols_to_use) > 1) {
        qrx <- qr(mat_x)
        pivot_idx <- qrx$pivot[seq_len(qrx$rank)]
        MatX <- mat_x[, pivot_idx, drop=FALSE]
        cols_used <- cols_to_use[pivot_idx]
      } else {
        MatX <- mat_x
        cols_used <- cols_to_use
      }
      
      BC_X <- cbind(BK_full_rank, MatX)
      fit_X <- lm.fit(BC_X, nu_hat)
      res_X <- nu_hat - as.vector(BC_X %*% fit_X$coefficients)
      Rn_X <- m * sum(res_X^2)
      
      p <- ncol(BC_X)
      lambda_BtB <- eigen(crossprod(BC_X), symmetric=TRUE, only.values=TRUE)$values
      ev_p <- (1 - (c_xmat^2) * (lambda_BtB^2))^2
      
      chi_part <- rchisq(500, df = K2 - p)
      Z_p <- matrix(rnorm(p * 500)^2, nrow = p, ncol = 500)
      B_X <- as.vector(crossprod(Z_p, ev_p)) + chi_part
      
      pval <- sum(B_X > as.numeric(Rn_X)) / 500
      return(list(pval = pval, MatX = MatX, cols = cols_used))
    }
    
    mod3 <- test_model_fit(tick3_col)
    mod2 <- test_model_fit(sort(c(tick3_col, tick2_col)))
    mod1 <- test_model_fit(sort(c(tick3_col, tick2_col, tick1_col)))
    
    final_mod <- NULL
    if(!is.null(mod3$MatX) && mod3$pval > 0.05) { final_mod <- mod3 } 
    else if (!is.null(mod2$MatX) && mod2$pval > 0.05) { final_mod <- mod2 } 
    else if (!is.null(mod1$MatX)) { final_mod <- mod1 }
    
    if (!is.null(final_mod) && !is.null(final_mod$MatX)) {
      keep_idx <- 1:length(final_mod$cols)
      
      while(length(keep_idx) > 0) {
        BC_X <- cbind(BK_full_rank, final_mod$MatX[, keep_idx, drop=FALSE])
        p_X <- ncol(BC_X)
        df_resid <- K2 - p_X
        
        if (df_resid <= 0) {
          worst_idx <- which.min(triad_scores[final_mod$cols[keep_idx]])
          keep_idx <- keep_idx[-worst_idx]
          next
        }
        
        fit_X <- lm.fit(BC_X, nu_hat)
        
        # Rank Deficiency Check
        triad_coefs <- fit_X$coefficients[K:p_X]
        if (any(is.na(triad_coefs))) {
          worst_idx <- which(is.na(triad_coefs))[1]
          keep_idx <- keep_idx[-worst_idx]
          next
        }
        
        sig2 <- sum(fit_X$residuals^2) / df_resid
        XtX_inv <- chol2inv(chol(crossprod(BC_X)))
        se <- sqrt(sig2 * diag(XtX_inv))
        
        triad_idx <- K:p_X
        t_stats <- fit_X$coefficients[triad_idx] / se[triad_idx]
        pvals <- 2 * pt(-abs(t_stats), df = df_resid)
        
        if (max(pvals) > 0.01) {
          keep_idx <- keep_idx[-which.max(pvals)]
        } else {
          break
        }
      }
      
      if (length(keep_idx) > 0) {
        final_cols_ftbs <- final_mod$cols[keep_idx]
        out_n_triad_ftbs <- length(final_cols_ftbs)
        out_true_model_select_ftbs <- check_span(final_cols_ftbs)
        
        final_MatX_ftbs <- final_mod$MatX[, keep_idx, drop=FALSE]
        BC_X_final_ftbs <- cbind(BK_full_rank, final_MatX_ftbs)
        final_hat_nu_ftbs <- lm.fit(BC_X_final_ftbs, nu_hat)$coefficients
      } else {
        out_n_triad_ftbs <- 0
        out_true_model_select_ftbs <- 0
        final_MatX_ftbs <- NULL
      }
    } else {
      out_n_triad_ftbs <- 0
      out_true_model_select_ftbs <- 0
      final_MatX_ftbs <- NULL
    }
    
    #=================================================#
    # DECISION-RELEVANT METRICS (Evaluated on Total Margin)
    #=================================================#
    
    # 1. True Model Total Prediction
    X_true <- cbind(BK_full_rank, CK_t)
    fit_true <- lm.fit(X_true, nu_hat)
    coef_true <- ifelse(is.na(fit_true$coefficients), 0, fit_true$coefficients)
    hat_nu_true_tot <- as.vector(X_true %*% coef_true)
    
    # 2. LASSO Total Prediction
    if (length(lasso_select) > 0) {
      hat_nu_lasso_tot <- hat_nu_linear + as.vector(as.matrix(CK[, lasso_select, drop=FALSE]) %*% lasso_coefs[lasso_select])
    } else {
      hat_nu_lasso_tot <- hat_nu_linear
    }
    
    # 3. Hybrid Total Prediction
    if (out_n_triad_hybrid > 0) {
      BC_X_hybrid <- cbind(BK_full_rank, final_MatX_hybrid[, keep_idx_hybrid, drop=FALSE])
      fit_hybrid <- lm.fit(BC_X_hybrid, nu_hat)
      coef_hybrid <- ifelse(is.na(fit_hybrid$coefficients), 0, fit_hybrid$coefficients)
      hat_nu_hybrid_tot <- as.vector(BC_X_hybrid %*% coef_hybrid)
    } else {
      hat_nu_hybrid_tot <- hat_nu_linear
    }
    
    # 4. Pure FTBS Total Prediction
    if(!is.null(final_MatX_ftbs)) {
      coef_ftbs <- ifelse(is.na(final_hat_nu_ftbs), 0, final_hat_nu_ftbs)
      hat_nu_ftbs_tot <- as.vector(BC_X_final_ftbs %*% coef_ftbs)
    } else {
      hat_nu_ftbs_tot <- hat_nu_linear
    }
    
    # 5. Full & 6. Reduced Total Predictions
    hat_nu_full_tot <- nu_hat 
    hat_nu_red_tot <- hat_nu_linear
    
    # ORIGINAL METRIC DEFINITION (SSE on the Total Margin)
    calc_metrics <- function(pred_tot) {
      c(
        sum((pred_tot - nu)^2),           
        sum(sign(round(pred_tot, 6)) != sign(round(nu, 6)))   
      )
    }
    
    met_true   <- calc_metrics(hat_nu_true_tot)
    met_lasso  <- calc_metrics(hat_nu_lasso_tot)
    met_ftbs   <- calc_metrics(hat_nu_ftbs_tot)
    met_hybrid <- calc_metrics(hat_nu_hybrid_tot)
    met_full   <- calc_metrics(hat_nu_full_tot)
    met_red    <- calc_metrics(hat_nu_red_tot)
    
    return(c(
      out_n_triad_lasso, out_true_model_select_lasso,
      out_n_triad_hybrid, out_true_model_select_hybrid,
      out_n_triad_ftbs, out_true_model_select_ftbs,
      met_true[1], met_lasso[1], met_ftbs[1], met_hybrid[1], met_full[1], met_red[1],
      met_true[2], met_lasso[2], met_ftbs[2], met_hybrid[2], met_full[2], met_red[2]
    ))
    
  }, error = function(err) {
    return(rep(NA, 18)) 
  })
  
}, mc.cores = num_cores)

#-----------------------------#
# 6. Aggregate Results & Output
#-----------------------------#
res_mat <- do.call(rbind, results_list)
colnames(res_mat) <- c(
  "n_triad_lasso", "true_mod_lasso", 
  "n_triad_hybrid", "true_mod_hybrid", 
  "n_triad_ftbs", "true_mod_ftbs", 
  "MSE_true", "MSE_lasso", "MSE_ftbs", "MSE_hybrid", "MSE_full", "MSE_red",
  "RE_true", "RE_lasso", "RE_ftbs", "RE_hybrid", "RE_full", "RE_red"
)

cat("\n==================================================\n")
cat("SIMULATION RESULTS ( K =", K, "| m =", m, "| Reps =", rep_sim, ")\n")
cat("==================================================\n")

cat("\n[TABLE 2]: True Model Recovery & Expected Size\n")
tab2 <- matrix(c(
  sum(res_mat[,"true_mod_lasso"] == n_true, na.rm=T) / rep_sim,  mean(res_mat[,"n_triad_lasso"], na.rm=T) / n_true,
  sum(res_mat[,"true_mod_hybrid"] == n_true, na.rm=T) / rep_sim, mean(res_mat[,"n_triad_hybrid"], na.rm=T) / n_true,
  sum(res_mat[,"true_mod_ftbs"] == n_true, na.rm=T) / rep_sim,   mean(res_mat[,"n_triad_ftbs"], na.rm=T) / n_true
), nrow = 3, byrow = T)
rownames(tab2) <- c("Standard LASSO", "Hybrid (LASSO + FTBS Pruning)", "Pure FTBS (FDR + Pruning)")
colnames(tab2) <- c("P(True Signal in Span)", paste0("E[Model Size / ", n_true, "]"))
print(round(tab2, 4))

cat("\n[TABLE 4]: Decision-Relevant Performance (Evaluated on Total Margin)\n")
tab4 <- data.frame(
  Method = c("True", "LASSO", "FTBS", "LASSO+FTBS", "Full", "Reduced"),
  SSE = sprintf("%.4f", c(  # Changed label from MSE_nu to SSE to match unscaled sum
    mean(res_mat[,"MSE_true"], na.rm=T),
    mean(res_mat[,"MSE_lasso"], na.rm=T),
    mean(res_mat[,"MSE_ftbs"], na.rm=T),
    mean(res_mat[,"MSE_hybrid"], na.rm=T),
    mean(res_mat[,"MSE_full"], na.rm=T),
    mean(res_mat[,"MSE_red"], na.rm=T)
  )),
  RE = sprintf("%.4f", c(
    mean(res_mat[,"RE_true"], na.rm=T),
    mean(res_mat[,"RE_lasso"], na.rm=T),
    mean(res_mat[,"RE_ftbs"], na.rm=T),
    mean(res_mat[,"RE_hybrid"], na.rm=T),
    mean(res_mat[,"RE_full"], na.rm=T),
    mean(res_mat[,"RE_red"], na.rm=T)
  )),
  Size_Ratio = c(
    "1.00",
    sprintf("%.2f", mean(res_mat[,"n_triad_lasso"], na.rm=T) / n_true),
    sprintf("%.2f", mean(res_mat[,"n_triad_ftbs"], na.rm=T) / n_true),
    sprintf("%.2f", mean(res_mat[,"n_triad_hybrid"], na.rm=T) / n_true),
    "-",
    "0.00"
  )
)
print(tab4, row.names = FALSE, right = TRUE)
cat("==================================================\n")