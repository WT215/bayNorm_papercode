load("E:/RNAseqProject/Illustrator_bayNorm/bayNorm_papercode/RealData/Tung_study/Load_Tung.RData")
#load DE functions
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

#bayNorm########
library(bayNorm)
bayNorm_mean_N1<-bayNorm(Data=N1_DAT,BETA_vec = efficiency[colnames(N1_DAT)],S=1000,mean_version = T)
bayNorm_mean_N2<-bayNorm(Data=N2_DAT,BETA_vec = efficiency[colnames(N2_DAT)],S=1000,mean_version = T)
bayNorm_mean_N3<-bayNorm(Data=N3_DAT,BETA_vec = efficiency[colnames(N3_DAT)],S=1000,mean_version = T)


library(bayNorm)
bayNorm_array_N1<-bayNorm(Data=N1_DAT,BETA_vec = efficiency[colnames(N1_DAT)],S=5)
bayNorm_array_N2<-bayNorm(Data=N2_DAT,BETA_vec = efficiency[colnames(N2_DAT)],S=5)
bayNorm_array_N3<-bayNorm(Data=N3_DAT,BETA_vec = efficiency[colnames(N3_DAT)],S=5)

#DE between individuals NA19101 and NA19239
library(abind)
qq<-abind(bayNorm_array_N2$Bay_array,bayNorm_array_N3$Bay_array,along=2)
dim(qq)
M_bay_mat<-SCnorm_runMAST3(Data=qq,NumCells = c(201,221))

#SAVER########
library(SAVER)
saver_N1<-saver(x=N1_DAT)
saver_N2<-saver(x=N2_DAT)
saver_N3<-saver(x=N3_DAT)

#generate 5 samples from SAVER posterior, in consistent with bayNorm
saver_N1_s<-sample.saver(saver_N1,rep=5)
saver_N2_s<-sample.saver(saver_N2,rep=5)
saver_N3_s<-sample.saver(saver_N3,rep=5)

library(abind)
saver_array_N1<-abind(saver_N1_s,along=3)
saver_array_N2<-abind(saver_N2_s,along=3)
saver_array_N3<-abind(saver_N3_s,along=3)
saver_array_N1[is.na(saver_array_N1)]<-0
saver_array_N2[is.na(saver_array_N2)]<-0
saver_array_N3[is.na(saver_array_N3)]<-0

#MAST_saver<-SCnorm_runMAST3(Data=cbind(saver_N2$estimate,saver_N3$estimate),NumCells = c(201,221))
qq<-abind(saver_array_N2,saver_array_N3,along=2)
#apply MAST on each one of 5 samples
MAST_saver_array_mat<-SCnorm_runMAST3(Data=qq,NumCells = c(201,221))


#####SCnorm#######
library(SCnorm)
DAT_ALL<-as.matrix(cbind(N1_DAT,N2_DAT,N3_DAT))
#set ditherCounts = TRUE since these are UMI data
scnorm_norm_ori<-SCnorm(Data=DAT_ALL,Conditions=c(rep(1,dim(N1_DAT)[2]),rep(2,dim(N2_DAT)[2]),rep(3,dim(N3_DAT)[2])),ditherCounts = TRUE)
scnorm_norm<-scnorm_norm_ori@metadata$NormalizedData
grs<-c(dim(N2_DAT)[2],dim(N3_DAT)[2])
M_scnorm = SCnorm_runMAST(Data=scnorm_norm[,colnames(cbind(N2_DAT,N3_DAT))], NumCells=as.numeric(grs))



##MAGIC#####
library(Rmagic)
load("C:/RAW_REAL/RAW_TUNG.RData")


MAGIC_N1<-run_magic(t(N1_DAT))
rownames(MAGIC_N1)<-rownames(t(N1_DAT))
colnames(MAGIC_N1)<-colnames(t(N1_DAT))
MAGIC_N1<-t(MAGIC_N1)

MAGIC_N2<-run_magic(t(N2_DAT))
rownames(MAGIC_N2)<-rownames(t(N2_DAT))
colnames(MAGIC_N2)<-colnames(t(N2_DAT))
MAGIC_N2<-t(MAGIC_N2)


MAGIC_N3<-run_magic(t(N3_DAT))
rownames(MAGIC_N3)<-rownames(t(N3_DAT))
colnames(MAGIC_N3)<-colnames(t(N3_DAT))
MAGIC_N3<-t(MAGIC_N3)

#MAGIC_TUNG<-cbind(MAGIC_N1,MAGIC_N2,MAGIC_N3)
MAGIC_TUNG <- run_magic(t(cbind(N1_DAT,N2_DAT,N3_DAT)))
rownames(MAGIC_TUNG)<-rownames(t(cbind(N1_DAT,N2_DAT,N3_DAT)))
colnames(MAGIC_TUNG)<-colnames(t(cbind(N1_DAT,N2_DAT,N3_DAT)))
MAGIC_TUNG<-t(MAGIC_TUNG)

#Run MAST on magic normalized data
M_MAGIC_tr<- SCnorm_runMAST(Data=cbind(MAGIC_N2,MAGIC_N3),NumCells=c(201,221))



#scImpute###########
library(scImpute)
scimpute_fun(Data=N1_DAT,Data_name = 'N1_DAT')
scimpute_fun(Data=N2_DAT,Data_name = 'N2_DAT')
scimpute_fun(Data=N3_DAT,Data_name = 'N3_DAT')

# #load scImpute
scImpute_N1_DAT<-readRDS(file="E:/RNAseqProject/tung2017batch/FINAL/scImpute/scImpute_N1_DAT.rds")
scImpute_N2_DAT<-readRDS(file="E:/RNAseqProject/tung2017batch/FINAL/scImpute/scImpute_N2_DAT.rds")
scImpute_N3_DAT<-readRDS(file="E:/RNAseqProject/tung2017batch/FINAL/scImpute/scImpute_N3_DAT.rds")

MAST_scimpute = SCnorm_runMAST(Data=cbind(scImpute_N2_DAT,scImpute_N3_DAT), NumCells=as.numeric(grs))



#Scaling method########
grs<-c(dim(N2_DAT)[2],dim(N3_DAT)[2])
RB_N1<-t(t(N1_DAT)/efficiency[colnames(N1_DAT)])
RB_N2<-t(t(N2_DAT)/efficiency[colnames(N2_DAT)])
RB_N3<-t(t(N3_DAT)/efficiency[colnames(N3_DAT)])

MAST_RB = SCnorm_runMAST(Data=cbind(RB_N2,RB_N3), NumCells=as.numeric(grs))

save.image(file='Tung_norms.RData')