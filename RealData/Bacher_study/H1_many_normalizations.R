load("RAW_INITIATE.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")

H1_data_comb<-cbind(H1_p24,H1_p96)[whichg_H1,]



#bayNorm##########
#Using spike-ins to estimate BETA
ERCC_BETA<-colSums(cbind(ERCC_H1_p24/20,ERCC_H1_p96/10))
#Assume that mean capture efficiency is 10%
ERCC_BETA<-ERCC_BETA/mean(ERCC_BETA)*0.1
Beta_H1<-ERCC_BETA

#Apply bayNorm on H1 datasets: mean of posterior
mBAY_H1<-bayNorm(Data=H1_data_comb,BETA_vec=Beta_H1,S=1000,parallel=T,NCores=5,FIX_MU = T,GR=F,Conditions=CONDITION_H1,BB_SIZE = T,mode_version = F,mean_version = T,UMI_sffl=c(20,10),Prior_type = 'GG',verbose = T)

#Apply bayNorm on H1 datasets: 20 samples
aBAY_H1<-bayNorm_sup(Data=H1_data_comb,Conditions=CONDITION_H1,UMI_sffl = c(20,10),S=20,PRIORS = mBAY_H1$PRIORS_LIST,BETA_vec = do.call(c,mBAY_H1$BETA))

library(abind)
qq<-abind(aBAY_H1$Bay_array_list,along=3)
MAST_bay_mat<-SCnorm_runMAST3(Data=qq,NumCells = c(92,92))


#SCnorm####
scnorm_out_H1<-SCnorm(Data=H1_data_comb,Conditions=CONDITION_H1)
M_SCnorm_H1<-SCnorm_runMAST(Data=scnorm_out_H1@metadata$NormalizedData,NumCells = c(92,92))

#MAGIC####
library(Rmagic)
rawinput<-as.matrix(cbind(H1_p24,H1_p96)[whichg_H1,])
MAGIC_H1 <- run_magic(t(rawinput))
rownames(MAGIC_H1)<-rownames(t(rawinput))
colnames(MAGIC_H1)<-colnames(t(rawinput))
MAGIC_H1<-t(MAGIC_H1)
M_magic_H1<-SCnorm_runMAST3(Data=MAGIC_H1,NumCells = c(92,92))


#Scaling####
ercc_totalh1<-colSums(cbind(ERCC_H1_p24,ERCC_H1_p96))
RB_norm_H1<-t(t(H1_data_comb)/(ercc_totalh1/mean(ercc_totalh1)*0.1))
M_RB_H1<-SCnorm_runMAST(Data=RB_norm_H1,NumCells = c(92,92))



save("H1_many_normalizations.RData")








