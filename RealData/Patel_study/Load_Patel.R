library(SummarizedExperiment)
#install_github("willtownes/patel2014gliohuman")
library(patel2014gliohuman) 
data(patel_tpm)
data(patel_counts) 

patel_pd = colData(patel_counts)
patel_counts = as.data.frame(as.matrix(assay(patel_counts)))
LABEL_instrument<-patel_pd$instrument
names(LABEL_instrument)<-colnames(patel_counts)


LABEL_instrument[is.na(LABEL_instrument)]


table(patel_pd$instrument)
dim(patel_counts)

patel_SC<-patel_counts[,patel_pd$sampleType=='SC']


Inputdat<-round(patel_SC/20)
QUANT<-quantile(colSums(Inputdat),c(0.1,0.9))


Inputdat<-Inputdat[,-which(colSums(Inputdat)<=QUANT[1] | colSums(Inputdat)>=QUANT[2])]

Inputdat<-Inputdat[-which(rowMeans(Inputdat)<1),]

library(bayNorm)
BETA_re<-BetaFun(Data=Inputdat,MeanBETA = 0.06)


BETA<-BETA_re$BETA
names(BETA)<-colnames(Inputdat)
summary(BETA)
length(which(BETA<0.01))

dim(Inputdat)
Inputdat<-Inputdat[,-which(BETA<0.01)]
BETA_use<-BETA[-which(BETA<0.01)]

all.equal(names(BETA_use),colnames(Inputdat))



library(bayNorm)
LABEL_INDIVIDUAL<-LABEL_instrument[colnames(Inputdat)]
table(LABEL_INDIVIDUAL)

bay_out<-bayNorm(Data=Inputdat,BETA_vec =BETA_use,mean_version = T,S=1000,BB_SIZE = T)

save.image('Patel2014_bay_out.RData')