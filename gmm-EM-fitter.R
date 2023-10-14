library(Rcpp)
library(RcppArmadillo)
library(data.table)

################################################################################
# EM algorithm on MULTIPLE samples
# Uses RcppArmadillo E-step
################################################################################
# Take command line argument
K <- commandArgs(trailingOnly=TRUE)
K <- strtoi(K)

################################################################################
# SOURCE function

sourceCpp("ompstepper.cpp")

################################################################################
# DATA

patients <- read.csv("P5samples100_01.csv")[,1]
setwd("/net/beegfs/scratch/prutten/H132_P5denovo")
#patients <- dir()
#patients <- patients[c(1:2,sample(c(3:195,200),80,replace=FALSE),196:199,201:214)]
M <- length(patients)

################################################################################
# Create subsample

Msub <- 50
#set.seed(1)
patientsubsample <- patients[sample(M,Msub,replace=FALSE)]
subsample <- matrix(NA, Msub * 1e4, 8)

cat("Model",K,"start building subsample from",Msub,"out of",M,"samples.","\n")

for (s in 1:Msub) {
  fcsFile <- patientsubsample[s]
  Y1 <- fread(file=fcsFile)
  Y1 <- as.matrix(Y1)
  
  a <- (s - 1) * 1e4 + 1
  b <- a + 1e4 - 1
  if (nrow(Y1) >= 1e4){
    subsample[a:b,] <- Y1[sample(nrow(Y1), 1e4, replace=FALSE),]
  } else {
    subsample[a:b,] <- Y1[sample(nrow(Y1), 1e4, replace=TRUE),]
  }
}

################################################################################
# Initialize parameters

p<-8
convergence_threshold <- 1e-3

means <- kmeans(subsample, K, iter.max=1000, algorithm="MacQueen")$centers
Sigmas <- t(matrix(rep(diag(p), K), p, K*p))
weights <- matrix(1/K, length(patients), K)

convergence <- 1
iter <- 1

################################################################################
# EM loop

while (convergence > convergence_threshold) {
  cat("Model",K,"iteration",iter,format(Sys.time(),usetz=TRUE),"\n")

  # Accumulators
  A0 <- matrix(0, length(patients), K) # rep(0,K)
  A1 <- matrix(0, K, p)
  A2 <- matrix(0, K*p, p)
  
  logL <- 0
  
  for (s in 1:M) {
    fcsFile <- patients[s]
    Y1 <- fread(file=fcsFile)
    Y1 <- as.matrix(Y1)
    #if (nrow(Y1) > 5e5) {
    #  set.seed(1)
    #  Y1 <- Y1[sample(1:nrow(Y1),5e5,replace=FALSE),]
    #}
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

################################################################################
# STORE RESULTS
setwd("/net/beegfs/scratch/prutten/results/run12")

write.csv(weights, paste("weights-H132-P5-100-samples-",toString(K),".csv",sep=''),row.names=FALSE)
write.csv(means, paste("means-H132-P5-100-samples-",toString(K),".csv",sep=''),row.names=FALSE)
write.csv(Sigmas, paste("covs-H132-P5-100-samples-",toString(K),".csv",sep=''),row.names=FALSE)




