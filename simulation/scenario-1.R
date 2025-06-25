library(MASS)
library(glmnet)
graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace
####################################
L_mat <- function(k,m) ## L matrix with equal observations each (i,j)
{
  Lm = matrix(0, nrow = choose(k,2)*m, ncol = choose(k,2))
  for (i in 1:choose(k,2)) {
    Lm[((i-1)*m+1):(i*m), i] = 1
  }
  return(Lm)
}
####################################
B_mat <- function(k) ## B matrix for complete graph
{
  M <- matrix(0, nrow = 1, ncol = k)
  for (b in 1:(k-1)) {
    bb = k-b
    M1 <- matrix(0, nrow = bb, ncol = k)
    ## 
    for (i in 1:bb) {
      M1[i,b] <- 1
      M1[i,i+b] <- -1
    }
    M <- rbind(M, M1)
  }
  return(M[-1,])
}
#############################################
### C matrix for the all distinct triads
C_mat<- function(k) 
{ 
  C1 = t(cbind(B_mat(k-1), diag(1, choose(k-1,2))))
  
  r0seq = seq(k-1, 3, by=-1)
  r1seq = 0
  for (j in 1:length(r0seq)) {
    r1seq[j] = sum(r0seq[1:j]) 
  }
  
  for (i in 1:(length(r0seq)-1)) { # i in seq(k-2, 3, by=-1)
    C1 = cbind(C1, rbind(matrix(0, nrow = r1seq[i], ncol = choose(r0seq[i]-1, 2)), 
                         t(cbind(B_mat(r0seq[i]-1), diag(1, choose(r0seq[i]-1, 2)))) ))
  }
  
  return(cbind(C1, c(rep(0, choose(k, 2)-3), 1,-1,1)))
}
#colSums(C_mat(10))
C_basis<- function(k) 
{ 
  return(t(cbind(B_mat(k-1), diag(1, choose(k-1,2)))))
}
#############################
## comparison pairs in nu vector components
########################
nu_vector <- function(K)
{
  nu_ij = matrix(0, nrow = 1, ncol = 2)
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      nu_ij = rbind(nu_ij, c(i,j))
    }
  }
  return(nu_ij[-1,])
}
#############################
## all cyclic triads list
########################
cyclic_triads<- function(K)
{
  cijk = matrix(0, nrow = 1, ncol = 3)
  for (i in 1:(K-2)) {
    for (j in (i+1):(K-1)) {
      for (k in (j+1):K) {
        cijk = rbind(cijk, c(i,j,k))
      }
    }
  }
  return(cijk[-1,])
}
############################################
## linear part parameter
mu_par <- function(K)
{
  nu1 = 0
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      nu1 = c(nu1, j-i)
    }
  }
  return(nu1[-1])
}
############################################
## simulation: scenario 2
############################################
## values of K and m
############################################
K=6; m=20
############################ the true parameter #############
# change coeffs and cyclic components for different nu 
a = 1; b = -1; c= 1; d =-1
nu = mu_par(K) + a*C_mat(K)[,1] + b*C_mat(K)[,8] # + c*C_mat(K)[,7] + d*C_mat(K)[,12]
############################################################
rep_sim=1000
n_triad = 0; true_model_select = 0; triad_tick= 0
n_triad_lasso=0; true_model_select_lasso=0
n_triad_fstep=0; true_model_select_fstep=0
n_triad_sts=0; true_model_select_sts=0

for (ai in 1:rep_sim) {
  
  ## ## data matrix
  bY = L_mat(K,m)%*%nu + rnorm(choose(K,2)*m, 0, 1)
  ########################################
  ## estimate of nu, nu_hat 
  ########################################
  nu_hat<- t(L_mat(K,m))%*%bY/m 
  #############################################################
  ## estimators of linear and cyclic by projecting nu_hat  ####
  #############################################################
  BK = B_mat(K)
  CK = C_mat(K)
  colnames(CK)= paste("v",1:ncol(CK), sep="")
  CKb = C_basis(K)
  hat_nu_linear = BK%*%ginv(t(BK)%*%BK)%*%t(BK)%*%nu_hat
  hat_nu_cyc = nu_hat-hat_nu_linear ## CK%*%ginv(t(CK)%*%CK)%*%t(CK)%*%nu_hat
  ######################
  ## variance of nu_hat
  ######################
  var_nu_hat = diag(choose(K,2), choose(K,2))
  var_hat_nu_cyc = CKb%*%solve(t(CKb)%*%CKb)%*%t(CKb)%*%var_nu_hat%*%
    t(CKb%*%solve(t(CKb)%*%CKb)%*%t(CKb))
  var_cyc_nu = diag(var_hat_nu_cyc)
  #######################
  ## testing v_ij \neq 0
  #######################
  nu_cyc_pval = 0
  
  for (i in 1:choose(K,2)) {
    nu_cyc_pval[i] = pnorm(abs(sqrt(m*choose(K,2))*hat_nu_cyc[i]/sqrt(var_cyc_nu[i])), 
                           mean=0, sd=1, lower.tail = F)
  }
  ## 
  nu_with_pval <- cbind(nu_vector(K), nu_cyc_pval)
  ## filtering significant nu_ij's using Bonferrani method
  hat_tr_nu <- nu_with_pval[nu_with_pval[,3]< 0.05/choose(K,2),]
  ## significant comparisons
  ######
  #hat_tr_nu[,1:2]
  # hat_tr_nu_order <- hat_tr_nu[order(hat_tr_nu$nu_cyc_pval, decreasing = F),]
  ## tick table
  all_triads = cyclic_triads(K)
  tick_table = matrix("", nrow = choose(K,2), ncol = choose(K,3))
  
  for (i in 1:choose(K,2)) {
    for (j in 1:choose(K,3)) {
      if(nu_with_pval[i,3]< 0.05/choose(K,2) & sum(nu_with_pval[i,1:2]%in%all_triads[j,])==2) 
      { tick_table[i,j] = paste("\u2713") }
    }
  }
  tick_table= as.data.frame(tick_table)
  colnames(tick_table) = apply(all_triads, 1, paste, collapse = ", ")
  rownames(tick_table) = apply(nu_vector(K), 1, paste, collapse = ", ")
  tick_table
  ## latex code of tick table
  library(xtable)
  xtable(tick_table)
  ###################################
  ## tick table: tick replaced with 1
  ###################################
  tick_table_10 = matrix(0, nrow = choose(K,2), ncol = choose(K,3))
  
  for (i in 1:choose(K,2)) {
    for (j in 1:choose(K,3)) {
      if(nu_with_pval[i,3]< 0.05/choose(K,2) & sum(nu_with_pval[i,1:2]%in%all_triads[j,])==2) 
      { tick_table_10[i,j] = 1 }
    }
  }
  tick_tab_colsum = colSums(tick_table_10)
  ## 3 tick columns
  tick3_col= which(tick_tab_colsum %in% 3)
  all_triads[which(tick_tab_colsum %in% 3),]
  ## 2 tick columns
  tick2_col= which(tick_tab_colsum %in% 2)
  all_triads[which(tick_tab_colsum %in% 2),]
  ## 1 tick columns
  tick1_col= which(tick_tab_colsum %in% 1)
  all_triads[which(tick_tab_colsum %in% 1),]
  
  ##############################################
  ### testing fit of model by our method ###
  ##########################################################################
  LB = L_mat(K,m)%*%BK ## LB
  L_Mat= L_mat(K,m)
  ############################################
  ## finding linearly independent cols 3
  ############################################
  mat3 = CK[,tick3_col]
  if(length(tick3_col)==1) 
  {Mat3= mat3} else{Mat3= mat3[, qr(mat3)$pivot[seq_len(qr(mat3)$rank)]]}
  
  LC3 = L_Mat%*%Mat3 ## LC
  ## xi matrix
  xmat = diag(1/choose(K,2), choose(K,2))
  
  ################
  ## estimate of nu with 3 tick columns only
  ##################
  hat_nu_003 = ginv(t(cbind(LB, LC3))%*%cbind(LB, LC3))%*%t(cbind(LB, LC3))%*%bY
  ## test statistic
  Rn_003 = t(nu_hat - cbind(BK,Mat3)%*%hat_nu_003)%*%
    diag(m, choose(K,2))%*%(nu_hat - cbind(BK,Mat3)%*%hat_nu_003)
  ## computing p-value intermediate quantities
  BC003 = cbind(BK,Mat3)
  xi003 = t(BC003)%*%xmat%*%BC003
  M003 = diag(1, choose(K,2))- BC003%*%xi003%*%t(BC003)%*%xmat
  psi_003 =sqrt(xmat)%*%M003%*%solve(xmat)%*%t(M003)%*%sqrt(xmat)
  #qr(psi003)$rank
  ev003 = eigen(psi_003)$values
  B003 = 0
  for (j in 1:500) {
    B003[j] = sum(ev003*(rnorm(choose(K,2), mean = 0, sd=1))^2)
  }
  ## estimated p-value  "high p-value shows good model fit"
  model_003_pval = length(B003[B003 > as.numeric(Rn_003)])/500
  ############################################
  ## finding linearly independent cols 3,2
  ############################################
  mat2 = CK[,sort(c(tick3_col,tick2_col))]
  if(length(c(tick3_col,tick2_col))==1)
  {Mat2=mat2} else{Mat2= mat2[, qr(mat2)$pivot[seq_len(qr(mat2)$rank)]]}
  LC2 = L_Mat%*%Mat2 ## LC
  ################
  ## estimate of nu with 3 and 2 ticks columns 
  ##################
  hat_nu_002 = ginv(t(cbind(LB, LC2))%*%cbind(LB, LC2))%*%t(cbind(LB, LC2))%*%bY
  ## test statistic
  Rn_002 = t(nu_hat - cbind(BK,Mat2)%*%hat_nu_002)%*%
    diag(m, choose(K,2))%*%(nu_hat - cbind(BK,Mat2)%*%hat_nu_002)
  ## computing p-value intermediate quantities
  BC002 = cbind(BK,Mat2)
  xi002 = t(BC002)%*%xmat%*%BC002
  M002 = diag(1, choose(K,2))- BC002%*%xi002%*%t(BC002)%*%xmat
  psi_002 =sqrt(xmat)%*%M002%*%solve(xmat)%*%t(M002)%*%sqrt(xmat)
  #qr(psi002)$rank
  ev002 = eigen(psi_002)$values
  B002 = 0
  for (j in 1:500) {
    B002[j] = sum(ev002*(rnorm(choose(K,2), mean = 0, sd=1))^2)
  }
  ## estimated p-value  "high p-value shows good model fit"
  model_002_pval = length(B002[B002 > as.numeric(Rn_002)])/500
  
  ############################################
  ## finding linearly independent cols 3,2,1
  ############################################
  mat1 = CK[,sort(c(tick3_col,tick2_col,tick1_col))]
  Mat1= mat1[, qr(mat1)$pivot[seq_len(qr(mat1)$rank)]]
  LC1 = L_Mat%*%Mat1 ## LC
  ################
  ## estimate of nu with 3,2 and 1 ticks columns 
  ##################
  hat_nu_001 = ginv(t(cbind(LB, LC1))%*%cbind(LB, LC1))%*%t(cbind(LB, LC1))%*%bY
  ## test statistic
  Rn_001 = t(nu_hat - cbind(BK,Mat1)%*%hat_nu_001)%*%
    diag(m, choose(K,2))%*%(nu_hat - cbind(BK,Mat1)%*%hat_nu_001)
  ## computing p-value intermediate quantities
  BC001 = cbind(BK,Mat1)
  xi001 = t(BC001)%*%xmat%*%BC001
  M001 = diag(1, choose(K,2))- BC001%*%xi001%*%t(BC001)%*%xmat
  psi_001 =sqrt(xmat)%*%M001%*%solve(xmat)%*%t(M001)%*%sqrt(xmat)
  #qr(psi001)$rank
  ev001 = eigen(psi_001)$values
  B001 = 0
  for (j in 1:500) {
    B001[j] = sum(ev001*(rnorm(choose(K,2), mean = 0, sd=1))^2)
  }
  ## estimated p-value  "high p-value shows good model fit"
  model_001_pval = length(B001[B001 > as.numeric(Rn_001)])/500
  ############################################
  ## outcome FTS #####
  ####################
  if(model_003_pval>0.05) 
  {n_triad[ai] = NCOL(Mat3); true_model_select[ai] = sum(tick3_col%in%1, tick3_col%in%8); triad_tick[ai] = 3} else if (model_002_pval>0.05)
  {n_triad[ai] = NCOL(Mat2); true_model_select[ai] = sum(c(tick3_col,tick2_col)%in%1, c(tick3_col,tick2_col)%in%8); triad_tick[ai] = 2} else
  {n_triad[ai] = NCOL(Mat1); true_model_select[ai] = sum(c(tick3_col,tick2_col,tick1_col)%in%1, c(tick3_col,tick2_col,tick1_col)%in%8); triad_tick[ai]= 1}
  ############################################
  ## outcome LASSO #####
  ######################
  LC = L_Mat%*%CK
  bY_adj = bY-L_Mat%*%hat_nu_linear
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet( LC, bY_adj, alpha = 1, intercept = F)
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  #view coefficients of best model
  best_model <- glmnet( LC, bY_adj, alpha = 1, lambda = best_lambda, intercept = F)
  lasso_select = which(coef(best_model)[-1]!=0)
  n_triad_lasso[ai]= length(lasso_select); 
  true_model_select_lasso[ai] = sum(lasso_select%in%1, lasso_select%in%8)
  ############################################
  ## outcome forward stepwise #####
  #################################
  bY_adj = bY-L_Mat%*%hat_nu_linear
  cov_1 =0
  for (i_1 in 1:NCOL(LC)) {
    cov_1[i_1] = summary(lm(bY_adj~ 0+ LC[,i_1]))$coefficients[,4]
  }
  cov_1_index= which(cov_1 == min(cov_1))
  #############
  i_max=1
  min_pval = 0
  while (min_pval < 0.05) {
    LC1 = LC[,-cov_1_index]
    cov_1 =0
    for (i_1 in 1:(NCOL(LC)-i_max)) {
      cov_1[i_1] = tail(summary(lm(bY_adj~ 0+ LC[,cov_1_index] + LC1[,i_1]))$coefficients[,4],1)
    }
    min_pval = min(cov_1)
    cov_0 = which(match(data.frame(LC), data.frame(LC1[,which(cov_1 == min(cov_1))]))==1)
    if(min_pval>0.05) {cov_2_index = cov_1_index} else {cov_2_index = c(cov_1_index, cov_0)}
    cov_1_index = cov_2_index
    #output
    i_max = i_max+1
  }
  fm.step_select = cov_1_index
  n_triad_fstep[ai]=length(fm.step_select); 
  true_model_select_fstep[ai]=sum(fm.step_select%in%1, fm.step_select%in%8)
  ##################################################################################
  ### stepwise triad selection
  ##################################################################################
  ## estimate of nu 
  #################
  bar_s_all<- t(L_mat(K,m))%*%bY/m 
  ###################################
  ## testing all the cyclic component (i,j,k)'s
  ##############################
  T_ijk = 0; p_ijk = 0
  
  for (i in 1:choose(K, 3)) {
    T_ijk[i] = sqrt(m/3)*t(C_mat(K)[,i])%*%bar_s_all ## test statistic
    ## p-value of the test
    p_ijk[i] = pnorm(abs(T_ijk[i]), mean=0, sd=1, lower.tail = F)
  }
  ##########################################################################
  ### testing fit of model by our method ###
  ##########################################################################
  pijk_0 = cbind(c(1:length(p_ijk)), p_ijk)
  pijk_order = pijk_0[order(pijk_0[,2], decreasing = F),]
  pijk_use = pijk_order[1:NROW(pijk_order),][,1]
  ##############
  LB = L_mat(K,m)%*%B_mat(K) ## LB
  LC = L_mat(K,m)%*%C_mat(K) ## LC
  B1 = B_mat(K); C1 = C_mat(K)
  xijk = diag(1/choose(K,2), choose(K,2))
  
  Rn_ijk = 0; model_ijk_pval = 0
  
  for (i in 1:length(pijk_use)) {
    ## estimate of nu if i^th column would be only true cyclic component 
    hat_nu_ijk = ginv(t(cbind(LB, LC[,pijk_use[1:i]]))%*%cbind(LB, LC[,pijk_use[1:i]]))%*%
      t(cbind(LB, LC[,pijk_use[1:i]]))%*%bY
    ## test statistic
    Rn_ijk[i] = t(bar_s_all - cbind(B1,C1[,pijk_use[1:i]])%*%hat_nu_ijk)%*%
      diag(m, choose(K,2))%*%(bar_s_all - cbind(B1,C1[,pijk_use[1:i]])%*%hat_nu_ijk)
    ## computing p-value intermediate quantities
    BCijk = cbind(B1,C1[,pijk_use[1:i]])
    xiijk = t(BCijk)%*%xijk%*%BCijk
    Mijk = diag(1, choose(K,2))- BCijk%*%xiijk%*%t(BCijk)%*%xijk
    psi_ijk =sqrt(xijk)%*%Mijk%*%solve(xijk)%*%t(Mijk)%*%sqrt(xijk)
    #qr(psiijk)$rank
    evijk = eigen(psi_ijk)$values
    Bijk = 0
    for (j in 1:500) {
      Bijk[j] = sum(evijk*(rnorm(choose(K,2), mean = 0, sd=1))^2)
    }
    ## estimated p-value
    model_ijk_pval[i] = length(Bijk[Bijk > as.numeric(Rn_ijk[i])])/500
    ## high p-value shows good model fit
  }
  ## output
  sts0 = which(model_ijk_pval > 0.05)[1]
  
  sts_select =pijk_use[1:sts0]
  n_triad_sts[ai]=length(sts_select) 
  true_model_select_sts[ai]=sum(sts_select%in%1, sts_select%in%8)
}

### output
fst_out= cbind(n_triad, true_model_select)
fst_out_true = fst_out[fst_out[,2]==2,] # true model selected part
nrow(fst_out_true)/rep_sim; mean(fst_out[,1])
## lasso
lasso_out = cbind(n_triad_lasso, true_model_select_lasso)
lasso_out_true = lasso_out[lasso_out[,2]==2,] # true model selected part
nrow(lasso_out_true)/rep_sim; mean(lasso_out[,1])
## forward stepwise
fstep_out = cbind(n_triad_fstep, true_model_select_fstep)
fstep_out_true = fstep_out[fstep_out[,2]==2,] # true model selected part
nrow(fstep_out_true)/rep_sim; mean(fstep_out[,1])
## sts
sts_out = cbind(n_triad_sts, true_model_select_sts)
sts_out_true = sts_out[sts_out[,2]==2,] # true model selected part
nrow(sts_out_true)/rep_sim; mean(sts_out[,1])
##
matrix(c(nrow(fstep_out_true)/rep_sim, mean(fstep_out[,1])/2,
         nrow(lasso_out_true)/rep_sim, mean(lasso_out[,1])/2,
         nrow(fst_out_true)/rep_sim, mean(fst_out[,1])/2,
         nrow(sts_out_true)/rep_sim, mean(sts_out[,1])/2), nrow = 4, byrow = T)

