source("E:/RNAseqProject/MANY_NORM_FUN.R")
load('Load_Torre.RData')

#bayNorm####
library(bayNorm)
bay_out<-bayNorm(Data=Torre_drop_sub,BETA_vec=BETA,S=5,NCores = 8)

#SAVER####
library(SAVER)
saver_out<-saver(x=Torre_drop_sub)
saver_s<-sample.saver(saver_out,rep=5)
library(abind)
saver_array<-abind(saver_s,along=3)

#MAGIC####
library(Rmagic)
MAGIC_Torre<-run_magic(data=t(Torre_drop_sub))
MAGIC_Torre<-t(MAGIC_Torre)
rownames(MAGIC_Torre)<-rownames(Torre_drop_sub)
colnames(MAGIC_Torre)<-colnames(Torre_drop_sub)
MAGIC_Torre<-as.matrix(MAGIC_Torre)


#scImpute####
scimpute_fun(Data=Torre_drop_sub,Data_name='scimpute_count')
scImpute_out<-readRDS("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/scImpute_outscimpute_count.rds")

#SCnorm####
library(SCnorm)
SCnorm_out<-SCnorm(Data=Torre_drop_sub,Conditions=rep(1,dim(Torre_drop_sub)[2]),ditherCounts = T)
SCnorm_dat<-SCnorm_out@metadata$NormalizedData

#Scaling####
RB_norm <- t(t(Torre_drop_sub) / bay_out$BETA)




save("Torre_many_normalizations.RData")



