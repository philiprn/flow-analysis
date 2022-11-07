library(flowCore)
library(flowAI)
library(PeacoQC)

################################################################################
# Pre-process raw FCS samples (and store)

setwd("/path/to/samples")
directory <- "/path/to/storage"

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
  ff <- RemoveDoublets(ff, channel1="FSC-A", channel2="FSC-H", nmad=4, verbose=FALSE, output="frame")
  
  # nrow(ff@exprs)
  
  ################################################################################
  c1<-which(colnames(ff@exprs)=="FITC-A")
  Y <- subset@exprs[,c1:(c1+7)]
  return(Y)
}

################################################################################

for (i in 1:length(dir())) {
  sam <- dir()[i]
  if (grep("P1",sam) == 1) {
    Y1 <- readFile(sam)
    fname <- paste(directory,"/",strtrim(sam,nchar(sam)-4),".csv",sep='')
    write.csv(Y1, fname, row.names=FALSE)
  }
}


