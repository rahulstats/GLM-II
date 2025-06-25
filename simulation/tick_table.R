library(MASS)
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
############################################
## parameter
nu_par <- function(K)
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
## values of K and m
############################################
K=6
############################ the true parameter #############
# change coeffs and cyclic components for different nu 
a = 1; b = -1; c= 1
nu = a*C_mat(K)[,8] + b*C_mat(K)[,16] + c*C_mat(K)[,18] ## scenario III
###########################################################
## 
nu_with <- cbind(nu_vector(K), nu)
## filtering significant nu_ij's using Bonferrani method
tr_nu <- nu_with[nu_with[,3]!=0,]
## tick table
############################################
all_triads = cyclic_triads(K)
tick_table = matrix("", nrow = choose(K,2), ncol = choose(K,3))

for (i in 1:choose(K,2)) {
  for (j in 1:choose(K,3)) {
    if(nu_with[i,3]!=0 & sum(nu_with[i,1:2]%in%all_triads[j,])==2) 
    { tick_table[i,j] = 123 #paste("\u2713") 
    }
  }
}
tick_table= as.data.frame(tick_table)
colnames(tick_table) = apply(all_triads, 1, paste, collapse = ", ")
rownames(tick_table) = apply(nu_vector(K), 1, paste, collapse = ", ")
#tick_table
## latex code of tick table
library(xtable)
xtable(tick_table)
