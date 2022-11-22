#paths
args <- commandArgs()
baseName <- args[6]

library(reticulate)
#use_python("C:/Users/PC/Desktop/mypython/venv/Scripts/python.exe", required = TRUE)
use_python("D:/Users/PC/anaconda3/envs/Objectdetection/python.exe", required = TRUE)
#C:/ProgramData/Microsoft/Windows/Start Menu/Programs/Anaconda3 (64-bit)
#source("runLeiden.R"); 
#source("scaleData.R"); 

infile <- paste(baseName, "H.txt", sep="/")
w10x.data <- read.table(file = infile,sep = '\t',row.names=1,header=T)
data.use <- w10x.data
#data.use <- scaleData(data.use, do.center = T)

SNN <- swne::CalcSNN(data.use, k = 20, prune.SNN = 1/15)

#out <- paste(baseName, "identity.txt", sep="/")
#write.table(as.matrix(idents),file = out, sep = '\t')

out <- paste(baseName, "identity.txt", sep="/")
write.table(as.matrix(SNN),file = out, sep = '\t')

