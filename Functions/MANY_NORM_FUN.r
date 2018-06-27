library(data.table)
#library(VennDiagram)
library(zoo)
library(statmod)
library(DESeq)
library(edgeR)

#scImpute
library(scImpute)
scimpute_fun<-function(Data,Data_name,labeled =FALSE,labels=NULL,...){
    saveRDS(Data,file=paste("~/",Data_name,".rds",sep=''))
    scimpute(paste("~/",Data_name,".rds",sep=''),infile='rds',outfile='rds',out_dir=paste("~/",Data_name,sep=''),ncores=1,labeled =labeled ,labels=labels,...)
    
}

## normalization methods
scran_norm <- function(x, sizes = c(20, 40, 60, 80, 100)) {
  library(scran)
  s <- computeSumFactors(x, sizes = sizes)
  counts <- t(t(x)/s*mean(s))
  n <- ncol(x)
  return(list(s=n*s/sum(s), counts=counts))
}

library(SCnorm)
scnorm_fun<-function(x,Conditions){
  library(SCnorm)
  DataNorm <- SCnorm(Data =x, Conditions = Conditions,PrintProgressPlots = F, FilterCellNum = 10, NCores=5)
}


rpm <- function(x) {
  s <- colSums(x)
  counts <- t(t(x)/colSums(x)*mean(s))
  n <- ncol(x)
  return(list(s=n*s/sum(s), counts=counts))
}
deseq <- function(x) {
  s <- estimateSizeFactorsForMatrix(x)
  counts <- t(t(x)/s*mean(s))
  n <- ncol(x)
  return(list(s=n*s/sum(s), counts=counts))
}
tmm <- function(x) {
  s <- edgeR::calcNormFactors(x,method= "TMM")*colSums(x)
  counts <- t(t(x)/s*mean(s))
  n <- ncol(x)
  return(list(s=n*s/sum(s), counts=counts))
}
div<-function(x,beta){
  counts <- t(t(x)/beta)
  return(list(s=beta, counts=counts))
}
