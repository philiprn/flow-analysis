library(Rcpp)
library(RcppArmadillo)
library(flowCore)

# function
sourceCpp("/Users/philip/Documents/postdoc/vumc/code/Rcpp/ares/getmeans.cpp")

# Sample file directory
setwd("/Users/philip/Documents/postdoc/vumc/data/MRD_NBM")

# Read a sample
fcsFile <- dir()[2]
subset <- read.FCS(fcsFile, transformation="linearize")
Y1 <- subset@exprs[,5:12]
Y1 <- as.matrix(Y1)

mu <- getm(Y1)

print(mu)






