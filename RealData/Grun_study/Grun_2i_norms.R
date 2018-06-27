load("Grun_2014_RAW.RData")
load('Grun_2014_RAW_serum.RData')
load("smFISH_norm_load.RData")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

CONDITION=c(rep('SC',dim(CountsUMI_SC)[2]),rep('P&S',dim(CountsUMI_RNA)[2]))
CONDITION_v2=c(rep(1,dim(CountsUMI_SC)[2]),rep(2,dim(CountsUMI_RNA)[2]))
names(BETA_ERCC)<-colnames(CountsUMI_2)



#use smFISH to estimate mean BETA######
FISH_mean<-unlist(lapply(Fish_2iplot_list,mean))[Fig3Gene]
SC_mean<-rowMeans(CountsUMI_SC[Fig3Gene,])[Fig3Gene]
RNA_mean<-rowMeans(CountsUMI_RNA[Fig3Gene,])[Fig3Gene]

#SC mean beta
lmm<-lm(SC_mean~FISH_mean)
lmm_summa<-summary(lmm)



#make Figure S11 (a)
plot(FISH_mean,SC_mean,pch=16,xlab='Mean expression of FISH counts',ylab='Mean expression of raw data',main='Data from Grün et al study (2i medium)',cex.main=0.95)
#text(FISH_mean,SC_mean,labels=names(FISH_mean),pos=)
legend('topleft',legend=c(paste('Estimated beta=',round(coef(lmm)[2],4)),paste('adj R2=',round(lmm_summa$adj.r.squared,4)),paste('cor=',round(cor(FISH_mean,SC_mean),4))),bty='n')
abline(lmm,lty=2)


summary(lmm)
MEANBETA_SC_2i<-coef(lmm)[2]
MEANBETA_SC_2i

#bayNorm####
library(bayNorm)
bayNorm_SC_2i<-bayNorm(Data=CountsUMI_SC,BETA_vec=BETA_ERCC_smFISH_SC_2i,Conditions=NULL,UMI_sffl = NULL,Prior_type = 'GG',mode_version = F,S=20,parallel=T,NCores=5,FIX_MU = T,GR=F,BB_SIZE = T,verbose=T)



#SAVER####
library(SAVER)
SAVER_SC_2i <- saver(CountsUMI_SC)
SAVER_SC_2i_samples<-sample.saver(SAVER_SC_2i, rep = 20, efficiency.known = FALSE, seed = NULL)
library(abind)
SAVER_SC_2i_array<-abind(SAVER_SC_2i_samples,along=3)

#scImpute####
library(scImpute)
scimpute_fun(Data=CountsUMI_SC,Data_name = 'scImpute_SC_2i')
scImpute_SC_2i<-readRDS("scImpute_SC_2i.rds")


#MAGIC####
MAGIC_2i <- run_magic(t(CountsUMI_SC))
MAGIC_2i<-as.matrix(MAGIC_2i)
rownames(MAGIC_2i)<-rownames(t(CountsUMI_SC))
colnames(MAGIC_2i)<-colnames(t(CountsUMI_SC))
MAGIC_2i<-t(MAGIC_2i)


#SCnorm####
library(SCnorm)
SCnorm_2i_ori<-SCnorm(x=CountsUMI_SC_2i,Conditions=c(rep(1,dim(CountsUMI_SC_2i)[2])),ditherCounts = T)


save.image('Grun_2i_norms.RData')