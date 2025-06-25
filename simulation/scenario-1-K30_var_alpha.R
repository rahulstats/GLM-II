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
## simulation: scenario 1
############################################
## values of K and m
############################################
K=30; m=100
alpha = 0.09
############################ the true parameter #############
# change coeffs and cyclic components for different nu 
a = 1; b = -1; c= 1; d =-1
nu = mu_par(K) + a*C_mat(K)[,1] + b*C_mat(K)[,56] 
#which(apply(cyclic_triads(30), 1, identical, c(1,4,5))==TRUE)
############################################################
rep_sim=1
n_triad = 0; true_model_select = 0; triad_tick= 0


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
  hat_tr_nu <- nu_with_pval[nu_with_pval[,3]< alpha/choose(K,2),]
  ## significant comparisons
  ######
  #hat_tr_nu[,1:2]
  # hat_tr_nu_order <- hat_tr_nu[order(hat_tr_nu$nu_cyc_pval, decreasing = F),]
  ## tick table
  all_triads = cyclic_triads(K)
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
  if(model_003_pval>alpha) 
  {n_triad[ai] = NCOL(Mat3); true_model_select[ai] = sum(tick3_col%in%16, tick3_col%in%50, tick3_col%in%66);
  triad_tick[ai] = 3} else if (model_002_pval>alpha)
  {n_triad[ai] = NCOL(Mat2); true_model_select[ai] = sum(c(tick3_col,tick2_col)%in%16, c(tick3_col,tick2_col)%in%50, 
  c(tick3_col,tick2_col)%in%66);triad_tick[ai] = 2} else
  {n_triad[ai] = NCOL(Mat1); true_model_select[ai] = 
   sum(c(tick3_col,tick2_col,tick1_col)%in%16, c(tick3_col,tick2_col,tick1_col)%in%50, c(tick3_col,tick2_col,tick1_col)%in%66);
  triad_tick[ai] = 1}
  
}
  
  