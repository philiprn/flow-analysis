library(Rcpp)
library(RcppArmadillo)
library(flowCore)
library(flowAI)
library(PeacoQC)

################################################################################
# Load GMM
setwd("/net/beegfs/scratch/prutten/results/run06/")

K<-30
p<-8

# Read model

weights <- read.csv("weights-H132-232-samples-30.csv")
means <- read.csv("means-H132-232-samples-30.csv")
Sigmas <- read.csv("covs-H132-232-samples-30.csv")

weights <- as.matrix(weights)
means <- as.matrix(means)
Sigmas <- as.matrix(Sigmas)

tot_weights <- colSums(weights)/dim(weights)[1]

################################################################################

sourceCpp("logclusterer.cpp")

################################################################################
# Function definition

readFile <- function(filename) {
  # Read file
  ff <- read.FCS(filename)
  
  # Remove "margin" events using FlowAI
  ff <- flow_auto_qc(ff,
                     remove_from='FM',
                     output = 1,
                     ChExcludeFM = 'Time',
                     ChExcludeFS = 'Time',
                     html_report=FALSE,
                     mini_report=FALSE,
                     folder_results=FALSE,
                     fcs_QC = FALSE)
  
  
  # Compensation
  ff <- compensate(ff, ff@description$SPILL)
  
  # Transformation
  ff <- transform(ff,transformList(colnames(ff@description$SPILL), 
                                   arcsinhTransform(a = 0, b = 1/150, c = 0)))
  
  # Remove anomalies
  ff <- flow_auto_qc(ff,
                     # deviationFR="MAD",
                     alphaFR=0.001,
                     remove_from='FS_FM',
                     output = 1,
                     # # ChExcludeFM = ignore_channels,
                     # ChExcludeFS = ignore_channels,
                     html_report=FALSE,
                     mini_report=FALSE,
                     folder_results=FALSE,
                     fcs_QC = FALSE)
  
  # nrow(ff@exprs)
  
  # Remove doublet cells using PeacoQC
  channel1 <- 'FSC-A'
  channel2 <- 'FSC-H'
  ff <- RemoveDoublets(ff, channel1="FSC-A", channel2="FSC-H", 
                       nmad=4, verbose=FALSE, output="frame")
  
  # nrow(ff@exprs)
  
  ##############################################################################
  c1<-which(colnames(ff@exprs)=="FITC-A")
  Y <- ff@exprs[,c1:(c1+7)]
  return(Y)
}

################################################################################

setwd("/net/beegfs/cfg/projects/mrd_mfc/H132_MRD")

tfiles<-read.csv("tfiles.csv")

N<-dim(tfiles)[1]

cat('Anlysing',N,"patients",'\n')

T0 <- matrix(NA, K, N)
T1 <- matrix(NA, K, N)
T2 <- matrix(NA, K, N)
T3 <- matrix(NA, K, N)

for (i in 1:N) {
  temp<-matrix(NA, K, 3)
  for (j in 1:length(dir())) {
    sam<-dir()[j]
    if (sam == tfiles[i,1]){
      ## Cluster sample #########################################
      Y1<-readFile(sam)
      dataClusters <- getClusters(Y1, tot_weights, means, Sigmas)
      for (k in 1:K) {
        T0[k,i] <- length(dataClusters[which(dataClusters == k)]) / length(dataClusters)
      }
      for (q in 1:3) {
        for (m in (j-25):(j+25)) {
          sam<-dir()[m]
          if (sam == tfiles[i,q+1]){
            ## Cluster sample #########################################
            Y1<-readFile(sam)
            dataClusters <- getClusters(Y1, tot_weights, means, Sigmas)
            for (k in 1:K) {
              temp[k,q] <- length(dataClusters[which(dataClusters == k)]) / length(dataClusters)
            }
            break
          }
        }
      }
    }
  }
  T1[,i] <- temp[,1]
  T2[,i] <- temp[,2]
  T3[,i] <- temp[,3]
}

################################################################################
# Store results

setwd("/net/beegfs/scratch/prutten/results/long01")

write.csv(T0, "denovo-clustering-4-tp.csv", row.names=FALSE)
write.csv(T1, "T1-clustering-4-tp.csv", row.names=FALSE)
write.csv(T2, "T2-clustering-4-tp.csv", row.names=FALSE)
write.csv(T3, "T3-clustering-4-tp.csv", row.names=FALSE)

