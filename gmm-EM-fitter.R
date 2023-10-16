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

sourceCpp("logstepper_logL_v2.cpp")

################################################################################
# DATA

setwd("/net/beegfs/scratch/prutten/H132_denovo")
#patients <- dir()
patients <- c(dir()[sample(1:211,105,replace=FALSE)],
	      dir()[sample(212:232,11,replace=FALSE)])

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
  Y1 <- fread(fcsFile)
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
    Y1 <- fread(fcsFile)
    Y1 <- as.matrix(Y1)
    if (nrow(Y1) > 5e5) {
      set.seed(1)
      Y1 <- Y1[sample(1:nrow(Y1),5e5,replace=FALSE),]
    }
    A0_ <- A0[s,]
    
    logL <- logL + logstep(Y1, weights[s,], means, Sigmas, A0_, A1, A2)
    
    A0[s,] <- A0_
    
    weights[s,] <- A0[s,] / nrow(Y1) 
  }
  
  if (iter > 1) {
    convergence <- abs((logL - logL_old) / logL_old)
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
setwd("/net/beegfs/scratch/prutten/results/run14")

write.csv(weights, paste("weights-H132-P1-232-samples-",toString(K),".csv",sep=''),row.names=FALSE)
write.csv(means, paste("means-H132-P1-232-samples-",toString(K),".csv",sep=''),row.names=FALSE)
write.csv(Sigmas, paste("covs-H132-P1-232-samples-",toString(K),".csv",sep=''),row.names=FALSE)

write.csv(patients,"trainset-H132-P1-232-samples.csv",row.names=FALSE)


