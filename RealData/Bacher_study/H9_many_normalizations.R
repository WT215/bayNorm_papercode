load("RAW_INITIATE.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")

H9_data_comb<-cbind(H9_p24,H9_p96)[whichg_H9,]



#bayNorm##########
#Using spike-ins to estimate BETA
ERCC_BETA<-colSums(cbind(ERCC_H9_p24/20,ERCC_H9_p96/10))
#Assume that mean capture efficiency is 10%
ERCC_BETA<-ERCC_BETA/mean(ERCC_BETA)*0.1
Beta_H9<-ERCC_BETA

#Apply bayNorm on H1 datasets: mean of posterior
mBAY_H9<-bayNorm(Data=H9_data_comb,BETA_vec=Beta_H9,S=1000,parallel=T,NCores=5,FIX_MU = T,GR=F,Conditions=CONDITION_H9,BB_SIZE = T,mode_version = F,mean_version = T,UMI_sffl=c(20,10),Prior_type = 'GG',verbose = T)

#Apply bayNorm on H1 datasets: 20 samples
aBAY_H9<-bayNorm_sup(Data=H9_data_comb,Conditions=CONDITION_H9,UMI_sffl = c(20,10),S=20,PRIORS = mBAY_H1$PRIORS_LIST,BETA_vec = do.call(c,mBAY_H9$BETA))

library(abind)
qq<-abind(aBAY_H9$Bay_array_list,along=3)
M_bay20_H9<-SCnorm_runMAST3(Data=qq,NumCells = c(91,91))


#SCnorm####
scnorm_out_H9<-SCnorm(Data=H9_data_comb,Conditions=CONDITION_H9)
M_SCnorm_H9<-SCnorm_runMAST(Data=scnorm_out_H9@metadata$NormalizedData,NumCells = c(91,91))

#MAGIC####
library(Rmagic)
rawinput<-as.matrix(cbind(H9_p24,H9_p96)[whichg_H9,])
MAGIC_H9 <- run_magic(t(rawinput))
rownames(MAGIC_H9)<-rownames(t(rawinput))
colnames(MAGIC_H9)<-colnames(t(rawinput))
MAGIC_H9<-t(MAGIC_H9)
M_magic_H9<-SCnorm_runMAST3(Data=MAGIC_H9,NumCells = c(91,91))


#Scaling####
ercc_totalh1<-colSums(cbind(ERCC_H9_p24,ERCC_H9_p96))
RB_norm_H9<-t(t(H1_data_comb)/(ercc_totalh1/mean(ercc_totalh1)*0.1))
M_RB_H9<-SCnorm_runMAST(Data=RB_norm_H9,NumCells = c(91,91))



save("H9_many_normalizations.RData")
