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

patients<-c()
for (i in 1:length(dir())) {
  sam<-dir()[i]
  if (length(grep("H132",sam)) != 1) {
    if (nchar(unlist(strsplit(sam," "))[2])==4){
      patients<-c(patients,unlist(strsplit(sam," "))[2])
    }
  }
}

patients<-unique(patients)
cat(length(patients),"patients",'\n')

cases<-0
for (i in 1:length(patients)) {
  denovo<-FALSE
  datesamples<-c()
  for (j in 1:length(dir())) {
    sam<-dir()[j]
    if (length(grep("H132",sam)) != 1) {
      if (nchar(unlist(strsplit(sam," "))[2])==4){
        if (unlist(strsplit(sam," "))[2] == patients[i]){
          if ((length(grep("BM",sam))==1) & (length(grep("P1",sam))==1)){
            if (length(grep("De Novo",sam)) == 1){
              denovo<-TRUE
            } else {
              datesamples<-c(datesamples,sam)
            }
          }
        } 
      }
    }
  }
  # Get number of cases
  if (length(datesamples) > 0){
    if ((length(datesamples) > 1) & denovo) {
      cases<-cases+1
    }
  }
}

cat(cases,"patients with 4 timepoints",'\n')
N <- cases

T0 <- matrix(NA, K, N)
T1 <- matrix(NA, K, N)
T2 <- matrix(NA, K, N)
T3 <- matrix(NA, K, N)

for (i in 1:length(patients)) {
  counter <- 1
  tempdates<-c()
  datesamples<-c()
  temp<-matrix(NA, K, 3)
  for (j in 1:dim(direc)[1]) {
    sam<-dir()[j]
    if (length(grep("H132",sam)) != 1) {
      if (nchar(unlist(strsplit(sam," "))[2])==4){
        if (unlist(strsplit(sam," "))[2] == patients[i]){
          if ((length(grep("BM",sam))==1) & (length(grep("P1",sam))==1)){
            if (length(grep("De Novo",sam)) == 1){
              ## Cluster sample #########################################
              Y1<-readFile(sam)
              dataClusters <- getClusters(Y1, tot_weights, means, Sigmas)
              for (k in 1:K) {
                T0[k,i] <- length(dataClusters[which(dataClusters == k)]) / length(dataClusters)
              }
            } else {
              datesamples<-c(datesamples,sam)
              tempdates<-c(tempdates,unlist(strsplit(sam," "))[3])
              ## Cluster sample ########################################
              Y1<-readFile(sam)
              dataClusters <- getClusters(Y1, tot_weights, means, Sigmas)
              for (k in 1:K) {
                temp[k,counter] <- length(dataClusters[which(dataClusters == k)]) / length(dataClusters)
              }
              counter <- counter + 1
            }
          }
        } 
      }
    }
  }
  # Sort samples by dates
  if (length(tempdates) > 0){
    tempdates<-as.Date(tempdates,format = "%d-%m-%Y")
    chrono<-order(tempdates)
    T1[,i] <- temp[,chrono[1]]
    T2[,i] <- temp[,chrono[2]]
    T3[,i] <- temp[,chrono[3]]
  }
}

################################################################################
# Store results

setwd("/net/beegfs/scratch/prutten/results/long01")

write.csv(T0, "denovo-clustering-4-tp.csv", row.names=FALSE)
write.csv(T1, "T1-clustering-4-tp.csv", row.names=FALSE)
write.csv(T2, "T2-clustering-4-tp.csv", row.names=FALSE)
write.csv(T3, "T3-clustering-4-tp.csv", row.names=FALSE)

