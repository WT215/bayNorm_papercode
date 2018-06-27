load('Load_Islam.RData')

library(bayNorm)
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")


#bayNorm####
BAY_LL_out<-bayNorm(Data=DAT_Jaakkola,BETA_vec=BETA,S=20,parallel=T,NCores=5,FIX_MU = T,GR=F,Conditions=CONDITION,BB_SIZE = T,mode_version = F,UMI_sffl=c(10,10),Prior_type = 'LL',verbose = T)



#MAGIC####
library(Rmagic)
MAGIC_islam <- run_magic(t(DAT_Jaakkola))
rownames(MAGIC_islam)<-rownames(t(DAT_Jaakkola))
colnames(MAGIC_islam)<-colnames(t(DAT_Jaakkola))
MAGIC_islam<-t(MAGIC_islam)
M_magic_islam<-SCnorm_runMAST(Data=MAGIC_islam,NumCells=c(48,44))

#SCnorm####
library(SCnorm)
SCnorm_out<-SCnorm(Data=DAT_Jaakkola,Conditions =CONDITION )
MAST_SCnorm<-SCnorm_runMAST(Data=SCnorm_out@metadata$NormalizedData,NumCells=c(48,44))

#scImpute#####
CONDITION<-c(rep(1,grs[1]),rep(2,grs[2]))
grs = table(substr(colnames(DAT_Jaakkola),1,2))
scimpute_fun(Data=DAT_Jaakkola,Data_name='scImpute_Jaakkolapars')

scImpute_Islam<-readRDS("scImpute_Jaakkolapars.rds")
M_scImpute_Islam<-SCnorm_runMAST(Data=scImpute_Islam,NumCells=grs)

#Scaling method####
sf_scran<-computeSumFactors(tttemp, sizes=c(20,30,40,50),positive=T)
BETA<-sf_scran/mean(sf_scran)*0.03
names(BETA)<-colnames(DAT_Jaakkola)
BETA[which(BETA==0)]=min(BETA[BETA>0])
sf_scran_rb<-computeSumFactors(DAT_Jaakkola, sizes=c(20,30,40,50),positive=T)
RB_norm<-t(t(DAT_Jaakkola))/(sf_scran_rb/mean(sf_scran_rb)*0.03)
MAST_RB<-SCnorm_runMAST(Data=RB_norm,NumCells=as.numeric(grs))

save.image('Islam_many_normalizations.RData')
