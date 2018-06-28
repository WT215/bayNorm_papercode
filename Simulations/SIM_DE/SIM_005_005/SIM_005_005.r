source("E:/RNAseqProject/QQsim_v2_SingleCellExperiment.R")
library(bayNorm)
load('Klein_bayNorm.RData')
Sim_List_Input$SCE<-NULL
Sim_List_Input$Est_params@counts.norm.TC<-matrix()




#######Begin simulation#######


Est_params<-Sim_List_Input$Est_params
Est_params@nGroups<-2
Est_params@groupCells=c(100,100)
Est_params@nGenes<-10000
#set mean BETA to be 5% in each group
Est_params@MeanBeta<-c(0.05,0.05)

#2000 out of 10000 genes were simulated to be DE in the first group
Est_params@de.prob<-c(0.2,0)
Est_params@de.downProb<-c(0.5,0)
Est_params@de.facLoc<-c(1,1)
Est_params@de.facScale<-c(0.5,0.5)

Est_params@out.prob
Est_params@mean.shape

Est_params@Beta<-numeric()
length(Est_params@Beta)


library(splatter)

SCE<-QQinitiate(Est_params)
SCE<-QQSimGeneMeans(SCE,Est_params)
#Group
SCE<-QQSimGroupDE(SCE,Est_params)
SCE<-QQSimGroupCellMeans(SCE, Est_params,Effect=F)

SCE<-QQSimBETA(SCE, Est_params)
SCE<-QQSimBCVMeans(SCE,Est_params)
SCE<-QQSimBinomial(SCE,Est_params)


dim(counts(SCE))
CONDITION<-SCE@colData@listData$GroupInd



DROP<-which(rowSums(counts(SCE))==0)

TRUE_LABEL<-rowData(SCE)$DEFacGroup1
TRUE_LABEL<-ifelse(TRUE_LABEL!=1,1,0)
TRUE_LABEL<-TRUE_LABEL[-DROP]
table(TRUE_LABEL)
#########True Parameters#####
TrueBeta1<-colData(SCE)$Beta[CONDITION==1]*Est_params@MeanBeta[1]
TrueBeta2<-colData(SCE)$Beta[CONDITION==2]*Est_params@MeanBeta[2]
summary(TrueBeta1)
summary(TrueBeta2)




TrueSIZES1<-1/(rowMeans(SCE@assays@.xData$data$BCV[,CONDITION==1])[-DROP])^2
TrueSIZES2<-1/(rowMeans(SCE@assays@.xData$data$BCV[,CONDITION==2])[-DROP])^2

TrueMU<-rowData(SCE)$BaseGeneMean[-DROP]



#begin bayNorm######
RAW_DAT<-counts(SCE)[-DROP,]
library(scran)
BETA_SCRAN_1<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1<-BETA_SCRAN_1/mean(BETA_SCRAN_1)*Est_params@MeanBeta[1]
BETA_SCRAN_2<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2<-BETA_SCRAN_2/mean(BETA_SCRAN_2)*Est_params@MeanBeta[2]

plot(BETA_SCRAN_2,TrueBeta2,pch=16)
abline(0,1)


library(bayNorm)
system.time(bayNorm_out<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,Prior_type = 'LL',S=10))


mbayNorm_out<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'LL')
gr<-table(CONDITION)
M_mbayNorm<-SCnorm_runMAST(Data=cbind(mbayNorm_out$Bay_mat_list$`Group 1`,mbayNorm_out$Bay_mat_list$`Group 2`),NumCells = gr)


library(abind)
qq<-abind(bayNorm_out$Bay_array_list,along=2)
dim(qq)
M_bay10<-SCnorm_runMAST3(Data=qq,NumCells=c(100,100))


load("E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005.RData")



######################other norms######

load("E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005.RData")



source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")


library(SCnorm)
scnorm_out_ori<-SCnorm(Data=RAW_DAT,Conditions=CONDITION,ditherCounts =T)
scnorm_out<-scnorm_out_ori@metadata$NormalizedData
M_scnorm<-SCnorm_runMAST(Data=scnorm_out,NumCells = gr)

#save(scnorm_out_ori,scnorm_out,M_scnorm,T_scnorm,file="E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005_SCnorm.RData")

#SAVER####
library(SAVER)
saver_temp1<-saver(x=RAW_DAT[,CONDITION==1])
saver_temp2<-saver(x=RAW_DAT[,CONDITION==2])
saver_out<-cbind(saver_temp1$estimate,saver_temp2$estimate)
M_saver<-SCnorm_runMAST(Data=saver_out,NumCells = c(100,100))


saver_s1<-sample.saver(saver_temp1,rep=10)
saver_s2<-sample.saver(saver_temp2,rep=10)
library(abind)
saver_array1<-abind(saver_s1,along=3)
sum(is.na(saver_array1))
saver_array2<-abind(saver_s2,along=3)
saver_array1[is.na(saver_array1)]<-0
saver_array2[is.na(saver_array2)]<-0

qq<-abind(saver_array1,saver_array2,along=2)
M_saver10<-SCnorm_runMAST3(Data=qq,NumCells = gr)

#save(saver_temp1,saver_temp2,saver_out,M_saver,M_saver10,saver_array1,saver_array2,file="E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005_saver_BETA.RData")


#scImpute####
library(scImpute)
scimpute_fun(Data=RAW_DAT,Data_name='DE_sim_005_005')
scImpute_out<-readRDS("E:/RNAseqProject/SIMULATION/SIM_005_005/DE_sim_005_005scimpute_count.rds")
gr<-table(CONDITION)
M_scImpute<-SCnorm_runMAST(Data=scImpute_out,NumCells = gr)


#Scaling method####
RB_out<-RAW_DAT/c(BETA_SCRAN_1,BETA_SCRAN_2)
M_RB<-SCnorm_runMAST(Data=RB_out,NumCells = gr)


save.image("E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005.RData")





###explore GG######
library(ROCR)
auc_fun<-function(pval,TRUELABEL)
{
    TRUE_LABEL_reverse<-TRUELABEL
    
    pred_MAST <- prediction((pval), TRUE_LABEL_reverse)
    perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
    
    auc_temp<-performance( pred_MAST, measure='auc' )
    auc_temp<-auc_temp@y.values[[1]]
    return(auc_temp)
}


load("E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005.RData")
library(bayNorm)
GG_bayNorm_out<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,Prior_type = 'GG',S=10)

GG_mbayNorm_out<-bayNorm_sup(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,mean_version = T,S=1000,PRIORS =GG_bayNorm_out$PRIORS_LIST )

library(abind)
qq<-abind(GG_bayNorm_out$Bay_array_list,along=2)
M_GG_bayNorm_out<-SCnorm_runMAST3(Data=qq,NumCells = c(100,100))
M_GG_mbayNorm_out<-SCnorm_runMAST3(Data=cbind(GG_mbayNorm_out$Bay_mat_list$`Group 1`,GG_mbayNorm_out$Bay_mat_list$`Group 2`),NumCells = c(100,100))

TRUE_LABEL_reverse<-TRUE_LABEL
TRUE_LABEL_reverse[TRUE_LABEL_reverse==1]=3
TRUE_LABEL_reverse[TRUE_LABEL_reverse==0]=1
TRUE_LABEL_reverse[TRUE_LABEL_reverse==3]=0


(auc_GG_bayNorm_out<-auc_fun(pval=apply(M_GG_bayNorm_out,1,median),TRUELABEL=TRUE_LABEL_reverse))
(auc_GG_mbayNorm_out<-auc_fun(pval=M_GG_mbayNorm_out$adjpval,TRUELABEL=TRUE_LABEL_reverse))

save(GG_bayNorm_out,GG_mbayNorm_out,M_GG_bayNorm_out,M_GG_mbayNorm_out,auc_GG_bayNorm_out,auc_GG_mbayNorm_out,file="E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_005_005.RData")



