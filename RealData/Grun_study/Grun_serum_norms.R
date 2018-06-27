load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/smFISH_norm_load.RData")
load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_RAW_serum.RData")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

CONDITION=c(rep('SC',dim(CountsUMI_SC)[2]),rep('P&S',dim(CountsUMI_RNA)[2]))
CONDITION_v2=c(rep(1,dim(CountsUMI_SC)[2]),rep(2,dim(CountsUMI_RNA)[2]))

#use smFISH to estimate mean BETA######
FISH_mean<-unlist(lapply(Fish_serumplot_list,mean))[Fig3Gene]
SC_mean<-rowMeans(CountsUMI_SC[Fig3Gene,])[Fig3Gene]
RNA_mean<-rowMeans(CountsUMI_RNA[Fig3Gene,])[Fig3Gene]

#SC mean beta
lmm<-lm(SC_mean~FISH_mean)
lmm_summa<-summary(lmm)



#generate figure for S11 (b)
pdf(file="E:/RNAseqProject/Illustrator_bayNorm/SUP/Grun_serumfit.pdf",width=4,height=4)
plot(FISH_mean,SC_mean,pch=16,xlab='Mean expression of FISH counts',ylab='Mean expression of raw data',main='Data from Grün et al study (serum medium)',cex.main=0.95)

legend('topleft',legend=c(paste('Estimated beta=',round(coef(lmm)[2],4)),paste('adj R2=',round(lmm_summa$adj.r.squared,4)),paste('cor=',round(cor(FISH_mean,SC_mean),4))),bty='n')
abline(lmm,lty=2)
dev.off()


summary(lmm)
MEANBETA_SC_serum<-coef(lmm)[2]

#RNA mean beta
lmm<-lm(RNA_mean~FISH_mean)
lmm_summa<-summary(lmm)
plot(FISH_mean,RNA_mean,pch=16,xlab='Mean expression of FISH counts',ylab='Mean expression of corresponding raw counts (Raw)',main=paste('Estimated beta=',round(coef(lmm)[2],4)))
legend('bottomright',legend=c(paste('adj R2=',round(lmm_summa$adj.r.squared,4)),paste('cor=',round(cor(FISH_mean,RNA_mean),4))))
abline(lmm)
summary(lmm)
MEANBETA_RNA_serum<-coef(lmm)[2]

BETA_ERCC_smFISH_SC_serum<-BETA_ERCC[CONDITION=='SC']/mean(BETA_ERCC[CONDITION=='SC'])*MEANBETA_SC_serum
summary(BETA_ERCC_smFISH_SC_serum)

BETA_ERCC_smFISH_RNA_serum<-BETA_ERCC[CONDITION=='P&S']/mean(BETA_ERCC[CONDITION=='P&S'])*MEANBETA_RNA_serum
summary(BETA_ERCC_smFISH_RNA_serum)


#bayNorm for serum####
library(bayNorm)
bayNorm_SC_serum<-bayNorm(Data=CountsUMI_SC,BETA_vec=BETA_ERCC_smFISH_SC_serum,Conditions=NULL,UMI_sffl = NULL,Prior_type = 'GG',mode_version = F,S=20,parallel=T,NCores=5,FIX_MU = T,GR=F,BB_SIZE = T,verbose=T)


#SAVER####
library(SAVER)
SAVER_SC_serum <- saver(CountsUMI_SC)
SAVER_SC_serum_samples<-sample.saver(SAVER_SC_serum, rep = 20, efficiency.known = FALSE, seed = NULL)
library(abind)
SAVER_SC_serum_array<-abind(SAVER_SC_serum_samples,along=3)

#MAGIC####
MAGIC_serum <- run_magic(t(CountsUMI_SC_serum))
MAGIC_serum<-as.matrix(MAGIC_serum)
rownames(MAGIC_serum)<-rownames(t(CountsUMI_SC_serum))
colnames(MAGIC_serum)<-colnames(t(CountsUMI_SC_serum))
MAGIC_serum<-t(MAGIC_serum)

#scImpute####
library(scImpute)
scimpute_fun(Data=CountsUMI_SC_serum,Data_name = 'scImpute_SC_serum')
scImpute_SC_serum<-readRDS("scImpute_SC_serum.rds")

#SCnorm####
library(SCnorm)
SCnorm_serum_ori<-SCnorm(x=CountsUMI_SC_serum,Conditions=c(rep(1,dim(CountsUMI_SC_serum)[2])),ditherCounts = T)

save.image('Grun_serum_norms.RData')
