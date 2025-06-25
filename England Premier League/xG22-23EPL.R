#Eng0=xgoals_pl_xlsx_Sheet2
#Eng1= as.matrix(Eng0)
#for (i in 1:19) {
 # Eng1[(19*i+1):(20*i),3]= (-1)*as.numeric(Eng1[(19*i+1):(20*i),3])
#}
#Eng1[,1:2]= t(apply(Eng1[,1:2], 1, sort))
#write.csv(Eng1, file = "Eng_xG_22-23.csv")
library(MASS)
library(glmnet)
graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace
setwd("~/Dropbox/UoH/g-lm2-simulations/updated/real data analysis")
Eng0=read.csv("xgoals.pl.xlsx - working.csv") ## EPL 22-23 xG
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
## nu vector components indices
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
######################
## real data analysis
#####################
K=20;m=2
bY=Eng0$xG
########################################
## estimate of nu, nu_hat 
########################################
nu_hat<- t(L_mat(K,m))%*%bY/m 
## estimating sigma^2
sigma2_hat= t(bY- L_mat(K,m)%*%nu_hat)%*%(bY- L_mat(K,m)%*%nu_hat)/length(bY)
############################################
## testing linear transitivity
###########################################
BK = B_mat(K)
## mu_hat
mu_hat = ginv(t(L_mat(K,m)%*%BK)%*% L_mat(K,m)%*%BK)%*%t(BK)%*% t(L_mat(K,m))%*%bY
##
D_m = diag(2, choose(K,2))
## U_n
U_lin = sqrt(D_m)%*%(diag(1, choose(K,2)) - BK%*% ginv(t(BK)%*%D_m%*%BK)%*%t(BK)%*%D_m )%*%nu_hat
## R_n
R_lin = t(U_lin)%*%U_lin/as.numeric(sigma2_hat)
## other matrices
xmat = diag(1/choose(K,2), choose(K,2))
xi00 = t(BK)%*%xmat%*%BK
M00 = diag(1, choose(K,2))- BK%*%ginv(xi00)%*%t(BK)%*%xmat
psi_00 = sqrt(xmat)%*%M00%*%solve(xmat)%*%t(M00)%*%sqrt(xmat)
#qr(psi00)$rank
ev00 = eigen(psi_00)$values
ev00 = round(ev00, 5)
## p-value for linearly transitive null
model_0_pval = pchisq(R_lin, 171, ncp = 0, lower.tail = F, log.p = FALSE)
model_0_pval ## pval
#############################################################
## estimators of linear and cyclic by projecting nu_hat  ####
#############################################################
CK = C_mat(K)
colnames(CK)= paste("v",1:ncol(CK), sep="")
CKb = C_basis(K)
hat_nu_linear = BK%*%ginv(t(BK)%*%BK)%*%t(BK)%*%nu_hat
hat_nu_cyc = nu_hat-hat_nu_linear ## CK%*%ginv(t(CK)%*%CK)%*%t(CK)%*%nu_hat
## norm of linear and cyclic parts
norm(nu_hat, type = "2")
norm(hat_nu_linear, type = "2")
norm(hat_nu_cyc, type = "2")
######################################################################
## verifying whether the cyclicality is purely by chance 
## assuming normal errors
#####################################################################
rep_sim =1000
nu = hat_nu_linear
norm_lin = 0; norm_cyc= 0
for (ai in 1:rep_sim) {
  ## ## data matrix
  bY1 = L_mat(K,m)%*%nu + rnorm(choose(K,2)*m, 0, sd=sigma2_hat)
  ## estimate of nu, nu_hat 
  nu_hat1<- t(L_mat(K,m))%*%bY1/m 
  ## estimators of linear and cyclic by projecting nu_hat  
  hat_nu_linear1 = BK%*%ginv(t(BK)%*%BK)%*%t(BK)%*%nu_hat1
  hat_nu_cyc1 = nu_hat1-hat_nu_linear1 ## CK%*%ginv(t(CK)%*%CK)%*%t(CK)%*%nu_hat
  ######################
  norm_lin[ai]= norm(hat_nu_linear1, type = "2")
  norm_cyc[ai]= norm(hat_nu_cyc1, type = "2")
}
##
mean(norm_lin); sd(norm_lin); mean(norm_lin < norm(hat_nu_linear, type = "2"))
mean(norm_cyc); sd(norm_cyc); mean(norm_cyc < norm(hat_nu_cyc, type = "2"))
########################################################################################
######################
## variance of nu_hat
######################
var_nu_hat = diag(choose(K,2), choose(K,2))
var_hat_nu_cyc = CKb%*%solve(t(CKb)%*%CKb)%*%t(CKb)%*%var_nu_hat%*%
  t(CKb%*%solve(t(CKb)%*%CKb)%*%t(CKb))
var_cyc_nu = as.numeric(sigma2_hat)*diag(var_hat_nu_cyc)
#######################
## testing v_ij \neq 0
#######################
nu_cyc_pval = 0
m=2
for (i in 1:choose(K,2)) {
  nu_cyc_pval[i] = pnorm(abs(sqrt(m*choose(K,2))*hat_nu_cyc[i]/sqrt(var_cyc_nu[i])), 
                         mean=0, sd=1, lower.tail = F)
}
## 
nu_with_pval <- cbind(nu_vector(K), nu_cyc_pval)
## filtering significant nu_ij's using Bonferrani method
hat_tr_nu <- nu_with_pval[nu_with_pval[,3]< 0.1/choose(K,2),]
## FDR
#hat_fdr_nu.p= p.adjust(nu_with_pval[,3], method = "fdr", n = length(nu_with_pval[,3]))
#hat_fdr_nu = cbind(nu_vector(K), hat_fdr_nu.p)
#hat_fdr.tr_nu <- hat_fdr_nu[hat_fdr_nu[,3]< 0.05,]
## significant comparisons
######
#hat_tr_nu[,1:2]
#hat_tr_nu_order <- hat_tr_nu[order(hat_tr_nu$nu_cyc_pval, decreasing = F),]
## tick table
all_triads = cyclic_triads(K)
tick_table = matrix("", nrow = choose(K,2), ncol = choose(K,3))

for (i in 1:choose(K,2)) {
  for (j in 1:choose(K,3)) {
    if(nu_with_pval[i,3]< 0.1/choose(K,2) & sum(nu_with_pval[i,1:2]%in%all_triads[j,])==2) 
    { tick_table[i,j] = paste("\u2713") }
  }
}
tick_table= as.data.frame(tick_table)
colnames(tick_table) = apply(all_triads, 1, paste, collapse = ", ")
rownames(tick_table) = apply(nu_vector(K), 1, paste, collapse = ", ")
#tick_table
## latex code of tick table
library(xtable)
xtable(tick_table)
###################################
## tick table: tick replaced with 1
###################################
tick_table_10 = matrix(0, nrow = choose(K,2), ncol = choose(K,3))

for (i in 1:choose(K,2)) {
  for (j in 1:choose(K,3)) {
    if(nu_with_pval[i,3]< 0.1/choose(K,2) & sum(nu_with_pval[i,1:2]%in%all_triads[j,])==2) 
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
## no triad with three tick
############################################
## xi matrix
xmat = diag(1/choose(K,2), choose(K,2))
############################################
## finding linearly independent cols 2
############################################
Mat2 = CK[,sort(c(tick2_col))]
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
psi_002 =as.numeric(sigma2_hat)*sqrt(xmat)%*%M002%*%solve(xmat)%*%t(M002)%*%sqrt(xmat)
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
psi_001 =as.numeric(sigma2_hat)*sqrt(xmat)%*%M001%*%solve(xmat)%*%t(M001)%*%sqrt(xmat)
#qr(psi001)$rank
ev001 = eigen(psi_001)$values
ev001 = round(ev001, 5)
B001 = 0
for (j in 1:500) {
  B001[j] = sum(ev001*(rnorm(choose(K,2), mean = 0, sd=1))^2)
}
## estimated p-value  "high p-value shows good model fit"
model_001_pval = length(B001[B001 > as.numeric(Rn_001)])/500
############################################
## outcome FTBS #####
####################
if (model_002_pval>0.05)
{n_triad  = NCOL(Mat2); triad_tick  = 2} else
{n_triad  = NCOL(Mat1); triad_tick = 1}
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
n_triad_lasso = length(lasso_select); n_triad_lasso
############################################
## outcome forward stepwise #####
#################################
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
n_triad_fstep =length(fm.step_select); n_triad_fstep
##################################################################################
### FSTS method
##################################################################################
## estimate of nu 
bar_s_all<- t(L_mat(K,m))%*%bY/m 
## testing all the cyclic component (i,j,k)'s
T_ijk = 0; p_ijk = 0

for (i in 1:choose(K, 3)) {
  T_ijk[i] = sqrt(m/3)*t(C_mat(K)[,i])%*%bar_s_all/sqrt(as.numeric(sigma2_hat)) ## test statistic
  ## p-value of the test
  p_ijk[i] = pnorm(abs(T_ijk[i]), mean=0, sd=1, lower.tail = F)
}
### testing fit  
pijk_0 = cbind(c(1:length(p_ijk)), p_ijk)
## plotting p-values lexicographically
#pdf(file="p_val_triads.pdf")
plot(pijk_0[,1], pijk_0[,2], pch = 20, type = "p", xlab = "Lexicographic order of the triads", ylab = "p-value",
     cex.lab=1.5, cex.axis=1.5 )
dev.off()

cdf1 = ecdf(pijk_0[,2])

plot(cdf1, cex=0)
## ordered p-val
pijk_order = pijk_0[order(pijk_0[,2], decreasing = F),]
pijk_use = pijk_order[1:NROW(pijk_order),][,1]
##############
LB = L_mat(K,m)%*%B_mat(K) ## LB
LC = L_mat(K,m)%*%C_mat(K) ## LC
B1 = B_mat(K); C1 = C_mat(K)
xijk = diag(1/choose(K,2), choose(K,2))

#Rn_ijk = 0; 
model_ijk_pval = 0; i=1
while (model_ijk_pval<= 0.05) {
  #for (i in 1:length(pijk_use)) {
  ## estimate of nu if i^th column would be only true cyclic component 
  hat_nu_ijk = ginv(t(cbind(LB, LC[,pijk_use[1:i]]))%*%cbind(LB, LC[,pijk_use[1:i]]))%*%
    t(cbind(LB, LC[,pijk_use[1:i]]))%*%bY
  ## test statistic
  Rn_ijk = t(bar_s_all - cbind(B1,C1[,pijk_use[1:i]])%*%hat_nu_ijk)%*%
    diag(m, choose(K,2))%*%(bar_s_all - cbind(B1,C1[,pijk_use[1:i]])%*%hat_nu_ijk)/as.numeric(sigma2_hat)
  ## computing p-value intermediate quantities
  BCijk = cbind(B1,C1[,pijk_use[1:i]])
  xiijk = t(BCijk)%*%xijk%*%BCijk
  Mijk = diag(1, choose(K,2))- BCijk%*%xiijk%*%t(BCijk)%*%xijk
  psi_ijk = sqrt(xijk)%*%Mijk%*%solve(xijk)%*%t(Mijk)%*%sqrt(xijk)
  #qr(psiijk)$rank
  evijk = eigen(psi_ijk)$values
  model_ijk_pval= pchisq(as.numeric(Rn_ijk), length(evijk), ncp = 0, lower.tail = F, log.p = FALSE)
  #Bijk = 0
  #for (j in 1:500) {
  #  Bijk[j] = sum(evijk*(rnorm(choose(K,2), mean = 0, sd=1))^2)
  #}
  ## estimated p-value
  #model_ijk_pval = length(Bijk[Bijk > as.numeric(Rn_ijk)])/500
  i=i+1
}
## output
sts0 = i-1 ## number of selected triads
sts_select =pijk_use[1:sts0]
n_triad_sts =sts0

all_triads[sts_select,]

############################################
## FTBS followed by FSR #####
#################################
cov_1 =0
LC0= LC[,sort(c(tick3_col,tick2_col,tick1_col))]
for (i_1 in 1:NCOL(LC0)) {
  cov_1[i_1] = summary(lm(bY_adj~ 0+ LC0[,i_1]))$coefficients[,4]
}
cov_1_index= which(cov_1 == min(cov_1))
#############
i_max=1
min_pval = 0
while (min_pval < 0.05) {
  LC1 = LC0[,-cov_1_index]
  cov_1 =0
  for (i_1 in 1:(NCOL(LC0)-i_max)) {
    cov_1[i_1] = tail(summary(lm(bY_adj~ 0+ LC0[,cov_1_index] + LC1[,i_1]))$coefficients[,4],1)
  }
  min_pval = min(cov_1)
  cov_0 = which(match(data.frame(LC0), data.frame(LC1[,which(cov_1 == min(cov_1))]))==1)
  if(min_pval>0.05) {cov_2_index = cov_1_index} else {cov_2_index = c(cov_1_index, cov_0)}
  cov_1_index = cov_2_index
  #output
  i_max = i_max+1
}
ftbs.fsr_select = cov_1_index
n_triad_ftbs.fsr =length(ftbs.fsr_select)

ftbs_fsr_selected_tr = 0
for (wh in 1:length(cov_1_index)) {
  ftbs_fsr_selected_tr[wh]= which(match(data.frame(LC) ,data.frame(LC0[,cov_1_index[wh]]))==1)
}

cbind(all_triads[ftbs_fsr_selected_tr,], pijk_0[ftbs_fsr_selected_tr,2])
### FTBS followed by FSR selected model
model_final=lm(bY~ 0+ LB+ LC[,ftbs_fsr_selected_tr])
coeff_vec= as.vector(model_final$coefficients)
coeff_vec[is.na(coeff_vec)] <- 0
##
final_nu_hat =  cbind(BK, CK[,ftbs_fsr_selected_tr])%*%coeff_vec

cbind(nu_vector(K), final_nu_hat)

final_nu_hat_linear = ginv(t(BK)%*%BK)%*%t(BK)%*%final_nu_hat

# hat_mu_s
hat_mu_s = ginv(t(LB)%*%LB)%*%t(LB)%*%(bY-L_mat(K,m)%*%CK[,ftbs_fsr_selected_tr]%*%coeff_vec[-c(1:K)])

#######################
## AIC BIC table
########################
## FSR
fsr_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,fm.step_select])
AIC(fsr_lm)
BIC(fsr_lm)
#FTBS
ftbs_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,sort(c(tick3_col,tick2_col,tick1_col))])
AIC(ftbs_lm)
BIC(ftbs_lm)
##FSTS
fsts_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,sts_select])
AIC(fsts_lm)
BIC(fsts_lm)
##FTBS+FSR
ftbs.fsr_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,ftbs_fsr_selected_tr])
AIC(ftbs.fsr_lm)
BIC(ftbs.fsr_lm)

matrix(c("FSR", length(fm.step_select), AIC(fsr_lm), BIC(fsr_lm),
         "FSTS", length(sts_select), AIC(fsts_lm), BIC(fsts_lm),
         "FTBS", length(c(tick3_col,tick2_col,tick1_col)), AIC(ftbs_lm), BIC(ftbs_lm),
         "FSTS-FSR", length(ftbs_fsr_selected_tr), AIC(ftbs.fsr_lm), BIC(ftbs.fsr_lm)), ncol=4, byrow = T)

##############################################################
############################# other analysis

initial_nu= cbind(nu_vector(K), nu_hat)

fitted_nu= cbind(nu_vector(K), final_nu_hat)

cbind(all_triads[pijk_order[,1],], pijk_order[,2])

ftbs_fsr_selected_tr = 0
for (wh in 1:length(cov_1_index)) {
  ftbs_fsr_selected_tr[wh]= which(match(data.frame(LC) ,data.frame(LC0[,cov_1_index[wh]]))==1)
}

cbind(all_triads[ftbs_fsr_selected_tr,], pijk_0[ftbs_fsr_selected_tr,2])

#################################################
## number of "lost" and "won" teams for each
#################################################
count_neg_nu=0

for (j in 1:20) {
  filter_nu=0
  for (i in 1:nrow(initial_nu)) {
    if(initial_nu[i,1]==j || initial_nu[i,2]==j) {filter_nu = c(filter_nu,initial_nu[i,3])}
  }
  
  filter_nu[1:j]= -1*filter_nu[1:j]
  filter_nu= filter_nu[-1]
  
  count_neg_nu[j]= sum(filter_nu<0)
}

loss_win_table= cbind(1:20, 19-count_neg_nu)

xtable(loss_win_table)

####################################
### rank set table.  ###
#####################################
Ai= list();Bi=list();Ei=list()

for (j in 1:20) {
  filter_nu=0
  for (i in 1:nrow(initial_nu)) {
    if(initial_nu[i,1]==j || initial_nu[i,2]==j) {filter_nu = c(filter_nu,initial_nu[i,3])}
  }
  
  filter_nu[1:j]= -1*filter_nu[1:j]
  filter_nu = filter_nu[-1]
  if(j==1) {filter_nu= c(0,filter_nu)} else 
  if(j==20) {filter_nu= c(filter_nu,0)} else  {filter_nu= c(filter_nu[1:(j-1)],0, filter_nu[j:19])}
  
  Ai[[j]]= which(filter_nu<0)
  Bi[[j]]= which(filter_nu>0)
  Ei[[j]]= which(filter_nu==0)
}





#################################################
## number of "lost" and "won" teams for each
#################################################
count_neg_nu=0

for (j in 1:20) {
  filter_nu=0
  for (i in 1:nrow(fitted_nu)) {
    if(fitted_nu[i,1]==j || fitted_nu[i,2]==j) {filter_nu = c(filter_nu,fitted_nu[i,3])}
  }
  
  filter_nu[1:j]= -1*filter_nu[1:j]
  filter_nu= filter_nu[-1]
  
  count_neg_nu[j]= sum(filter_nu<0)
}

loss_win_table= cbind(1:20, count_neg_nu, 19-count_neg_nu)

xtable(loss_win_table)
#####################################################################################
## checking whether after removing items in significant triads model becomes transitive
###################################################################################
teams = unique(Eng0[,1])
cyclic_teams = unique(teams[all_triads[ftbs_fsr_selected_tr,]])
Eng01 = subset(Eng0, !Home %in% cyclic_teams)
Eng02 = subset(Eng01, !Away %in% cyclic_teams)
##
K2=K-length(cyclic_teams);m=2
bY02=Eng02$xG
########################################
## estimate of nu, nu_hat 
########################################
nu_hat02<- t(L_mat(K2,m))%*%bY02/m 
## estimating sigma^2
sigma2_hat02= t(bY02- L_mat(K2,m)%*%nu_hat02)%*%(bY02- L_mat(K2,m)%*%nu_hat02)/length(bY02)
############################################
## testing linear transitivity
###########################################
BK2 = B_mat(K2)
D_m2 = diag(2, choose(K2,2))
## U_n
U_lin2 = sqrt(D_m2)%*%(diag(1, choose(K2,2)) - BK2%*% ginv(t(BK2)%*%D_m2%*%BK2)%*%t(BK2)%*%D_m2 )%*%nu_hat02
## R_n
R_lin2 = t(U_lin2)%*%U_lin2/as.numeric(sigma2_hat02)
## other matrices
xmat = diag(1/choose(K2,2), choose(K2,2))
xi00 = t(BK2)%*%xmat%*%BK2
M00 = diag(1, choose(K2,2))- BK2%*%ginv(xi00)%*%t(BK2)%*%xmat
psi_002 = sqrt(xmat)%*%M00%*%solve(xmat)%*%t(M00)%*%sqrt(xmat)
#qr(psi00)$rank
ev002 = eigen(psi_002)$values
ev002 = round(ev002, 5)
## p-value for linearly transitive null
model_02_pval = pchisq(R_lin2, choose(K2-1,2), ncp = 0, lower.tail = F, log.p = FALSE)
model_02_pval ## pval
## removed cyclic items 
removed_nu= cbind(nu_vector(K2), nu_hat02)
## 
#################################################
## number of "won" teams for each
#################################################
count_neg_nu=0

for (j in 1:17) {
  filter_nu=0
  for (i in 1:nrow(nu_hat02)) {
    if(removed_nu[i,1]==j || removed_nu[i,2]==j) {filter_nu = c(filter_nu,removed_nu[i,3])}
  }
  
  filter_nu[1:j]= -1*filter_nu[1:j]
  filter_nu= filter_nu[-1]
  
  count_neg_nu[j]= sum(filter_nu<0)
}

loss_win_table= cbind(1:17, 17-count_neg_nu)
loss_win_table
## cyclic triad in nu_hat02
cbind(nu_vector(K2), nu_hat02)





##################################################################################
### FTBS followed by FSTS 
##################################################################################
ftbs_all_triads = all_triads[sort(c(tick3_col,tick2_col,tick1_col)),]
C_matK2 = C_mat(K)[,sort(c(tick3_col,tick2_col,tick1_col))]
## estimate of nu 
bar_s_all<- t(L_mat(K,m))%*%bY/m 
## testing all the cyclic component (i,j,k)'s
T_ijk2 = 0; p_ijk2 = 0

for (i in 1:NCOL(C_matK2)) {
  T_ijk2[i] = sqrt(m/3)*t(C_matK2[,i])%*%bar_s_all/sqrt(as.numeric(sigma2_hat)) ## test statistic
  ## p-value of the test
  p_ijk2[i] = pnorm(abs(T_ijk2[i]), mean=0, sd=1, lower.tail = F)
}
### testing fit  
pijk2_0 = cbind(c(1:length(p_ijk2)), p_ijk2)
## plotting p-values lexicographically
#pdf(file="p_val_triads.pdf")
plot(pijk2_0[,1], pijk2_0[,2], pch = 20, type = "p", xlab = "Lexicographic order of the triads", ylab = "p-value",
     cex.lab=1.5, cex.axis=1.5 )
#dev.off()

cdf2 = ecdf(pijk2_0[,2])

plot(cdf2, cex=0)
## ordered p-val
pijk2_order = pijk2_0[order(pijk2_0[,2], decreasing = F),]
pijk2_use = pijk2_order[1:NROW(pijk2_order),][,1]
##############
LB = L_mat(K,m)%*%B_mat(K) ## LB
LC2 = L_mat(K,m)%*%C_matK2 ## LC
B1 = B_mat(K); C12 = C_matK2
xijk = diag(1/choose(K,2), choose(K,2))

#Rn_ijk = 0; 
model_ijk_pval = 0; i=1
while (model_ijk_pval<= 0.05) {
  ## estimate of nu if i^th column would be only true cyclic component 
  hat_nu_ijk = ginv(t(cbind(LB, LC2[,pijk2_use[1:i]]))%*%cbind(LB, LC2[,pijk2_use[1:i]]))%*%
    t(cbind(LB, LC2[,pijk2_use[1:i]]))%*%bY
  ## test statistic
  Rn_ijk = t(bar_s_all - cbind(B1,C12[,pijk2_use[1:i]])%*%hat_nu_ijk)%*%
    diag(m, choose(K,2))%*%(bar_s_all - cbind(B1,C12[,pijk2_use[1:i]])%*%hat_nu_ijk)/as.numeric(sigma2_hat)
  ## computing p-value intermediate quantities
  BCijk = cbind(B1,C12[,pijk2_use[1:i]])
  xiijk = t(BCijk)%*%xijk%*%BCijk
  Mijk = diag(1, choose(K,2))- BCijk%*%xiijk%*%t(BCijk)%*%xijk
  psi_ijk = sqrt(xijk)%*%Mijk%*%solve(xijk)%*%t(Mijk)%*%sqrt(xijk)
  #qr(psiijk)$rank
  evijk = eigen(psi_ijk)$values
  model_ijk_pval= pchisq(as.numeric(Rn_ijk), length(evijk), ncp = 0, lower.tail = F, log.p = FALSE)
  i=i+1
}
## output
n_triad_ftbs_sts = i-1 ## number of selected triads
ftbs_sts_select =pijk_use[1:n_triad_ftbs_sts ]
##
ftbs_sts_triads= all_triads[ftbs_sts_select,]
##
ftbs_fsts_selected_tr = 0
for (wh in 1:n_triad_ftbs_sts) {
  ftbs_fsts_selected_tr[wh]= which(apply(all_triads, 1, identical, ftbs_sts_triads[wh,]))
}
##
cbind(all_triads[ftbs_fsts_selected_tr,], pijk_0[ftbs_fsts_selected_tr,2])
##
### FTBS followed by FSTS selected model
model_final=lm(bY~ 0+ LB+ LC[,ftbs_fsts_selected_tr])
coeff_vec= as.vector(model_final$coefficients)
coeff_vec[is.na(coeff_vec)] <- 0
##
final_nu_hat =  cbind(BK, CK[,ftbs_fsts_selected_tr])%*%coeff_vec

cbind(nu_vector(K), final_nu_hat)

final_nu_hat_linear = ginv(t(BK)%*%BK)%*%t(BK)%*%final_nu_hat

# hat_mu_s
hat_mu_s = ginv(t(LB)%*%LB)%*%t(LB)%*%(bY-L_mat(K,m)%*%CK[,ftbs_fsts_selected_tr]%*%coeff_vec[-c(1:K)])
##
#####################################################################################
## checking whether after removing items in significant triads model becomes transitive
###################################################################################
teams = unique(Eng0[,1])
cyclic_teams = unique(teams[all_triads[ftbs_fsts_selected_tr,]])
Eng01 = subset(Eng0, !Home %in% cyclic_teams)
Eng02 = subset(Eng01, !Away %in% cyclic_teams)
##
K2=K-length(cyclic_teams);m=2
bY02=Eng02$xG
########################################
## estimate of nu, nu_hat 
########################################
nu_hat02<- t(L_mat(K2,m))%*%bY02/m 
## estimating sigma^2
sigma2_hat02= t(bY02- L_mat(K2,m)%*%nu_hat02)%*%(bY02- L_mat(K2,m)%*%nu_hat02)/length(bY02)
############################################
## testing linear transitivity
###########################################
BK2 = B_mat(K2)
D_m2 = diag(2, choose(K2,2))
## U_n
U_lin2 = sqrt(D_m2)%*%(diag(1, choose(K2,2)) - BK2%*% ginv(t(BK2)%*%D_m2%*%BK2)%*%t(BK2)%*%D_m2 )%*%nu_hat02
## R_n
R_lin2 = t(U_lin2)%*%U_lin2/as.numeric(sigma2_hat02)
## other matrices
xmat = diag(1/choose(K2,2), choose(K2,2))
xi00 = t(BK2)%*%xmat%*%BK2
M00 = diag(1, choose(K2,2))- BK2%*%ginv(xi00)%*%t(BK2)%*%xmat
psi_002 = sqrt(xmat)%*%M00%*%solve(xmat)%*%t(M00)%*%sqrt(xmat)
#qr(psi00)$rank
ev002 = eigen(psi_002)$values
ev002 = round(ev002, 5)
## p-value for linearly transitive null
model_02_pval = pchisq(R_lin2, choose(K2-1,2), ncp = 0, lower.tail = F, log.p = FALSE)
model_02_pval ## pval
## removed cyclic items 
removed_nu= cbind(nu_vector(K2), nu_hat02)
## 

#######################
## AIC BIC table
########################
## FSR
fsr_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,fm.step_select])
AIC(fsr_lm)
BIC(fsr_lm)
#FTBS
ftbs_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,sort(c(tick3_col,tick2_col,tick1_col))])
AIC(ftbs_lm)
BIC(ftbs_lm)
##FSTS
fsts_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,sts_select])
AIC(fsts_lm)
BIC(fsts_lm)
##FTBS+FSR
ftbs.fsr_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,ftbs_fsr_selected_tr])
AIC(ftbs.fsr_lm)
BIC(ftbs.fsr_lm)

##FTBS+FSTS
ftbs.fsts_lm= lm(bY~ 0+ LB+ (L_Mat%*%CK)[,ftbs_fsts_selected_tr])
AIC(ftbs.fsts_lm)
BIC(ftbs.fsts_lm)


matrix(c("FSR", length(fm.step_select), AIC(fsr_lm), BIC(fsr_lm),
         "FSTS", length(sts_select), AIC(fsts_lm), BIC(fsts_lm),
         "FTBS", length(c(tick3_col,tick2_col,tick1_col)), AIC(ftbs_lm), BIC(ftbs_lm),
         "FTBS-FSR", length(ftbs_fsr_selected_tr), AIC(ftbs.fsr_lm), BIC(ftbs.fsr_lm),
         "FTBS-FSTS", length(ftbs_fsts_selected_tr),AIC(ftbs.fsts_lm), BIC(ftbs.fsts_lm)), ncol=4, byrow = T)
