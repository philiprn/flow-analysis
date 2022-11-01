library(Rcpp)
library(RcppArmadillo)
library(flowCore)

################################################################################
# EM algorithm on MULTIPLE files
# Uses RcppArmadillo E-step
################################################################################
# Take command line argument
K <- commandArgs(trailingOnly=TRUE)
K <- strtoi(K)

################################################################################
# DATA

directory <- "/net/beegfs/cfg/projects/mrd_mfc/scratch-share/AML_NBM/patients/"
patients <- dir(directory)
M <- length(patients)

################################################################################
# SOURCE function

sourceCpp("/Users/philip/Documents/postdoc/vumc/code/Rcpp/logstepper_logL.cpp")

################################################################################
# Initialize parameters

p<-8
convergence_threshold <- 5e-3

subsample <- matrix(NA, 450000, 8)

for (s in 1:M) {
  fcsFile <- dir(paste(directory,patients[s],sep=''))[2]
  subset <- read.FCS(fcsFile, transformation="linearize")
  c1<-which(colnames(subset@exprs)=="FITC-A")
  Y1 <- subset@exprs[,c1:(c1+7)]
  Y1 <- as.matrix(Y1)
  
  a <- (s - 1) * 10000 + 1
  b <- a + 10000 - 1
  subsample[a:b,] <- Y1[sample(nrow(Y1), 10000, replace=FALSE),]
}

means <- kmeans(subsample, K, iter.max=1000, nstart=10, algorithm="MacQueen")$centers
Sigmas <- t(matrix(rep(diag(p), K), p, K*p))
weights <- matrix(1/K, length(patients), K)

convergence <- 1
iter <- 1

while (convergence > convergence_threshold) {
  # Accumulators
  A0 <- matrix(0, length(patients), K) # rep(0,K)
  A1 <- matrix(0, K, p)
  A2 <- matrix(0, K*p, p)
  
  logL <- 0
  
  for (s in 1:M) {
    fcsFile <- dir(paste(directory,patients[s],sep=''))[2]
    subset <- read.FCS(fcsFile, transformation="linearize")
    c1<-which(colnames(subset@exprs)=="FITC-A")
    Y1 <- subset@exprs[,c1:(c1+7)]
    Y1 <- as.matrix(Y1)
    
    A0_ <- A0[s,]
    
    logL <- logL + logstep(Y1, weights[s,], means, Sigmas, A0_, A1, A2)
    
    A0[s,] <- A0_
    
    weights[s,] <- A0[s,] / nrow(Y1) 
  }
  
  if (iter > 1) {
    convergence <- (logL_old - logL) / logL_old
  }
  
  logL_old <- logL
  iter <- iter + 1
  
  # M-step
  for (k in 1:K) {
    means[k,] <- A1[k,] / colSums(A0)[k]
    Sigmas[(k-1)*p + c(1:p),] <- A2[(k-1)*p + c(1:p),] / colSums(A0)[k] - t(tcrossprod(means[k,], means[k,]))
  }
}

# STORE RESULTS
setwd("~/results/")

write.csv(weights, paste("weights-MRD_NBM-48-samples-",toString(K),".csv",sep=''))
write.csv(means, paste("means-MRD_NBM-48-samples-",toString(K),".csv",sep=''))
write.csv(Sigmas, paste("covs-MRD_NBM-48-samples-",toString(K),".csv",sep=''))




