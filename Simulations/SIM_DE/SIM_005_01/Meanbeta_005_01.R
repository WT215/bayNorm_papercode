load("E:/RNAseqProject/SIMULATION/SIM_005_01/SIM_005_01.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")
library(bayNorm)

library(scran)
BETA_SCRAN_1_0025<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1_0025<-BETA_SCRAN_1_0025/mean(BETA_SCRAN_1_0025)*0.025
BETA_SCRAN_2_005<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2_005<-BETA_SCRAN_2_005/mean(BETA_SCRAN_2_005)*0.05



##### 0.025  0.05#########

mbayNorm_out_0025005<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1_0025,BETA_SCRAN_2_005),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'LL')

gr<-table(CONDITION)
M_mbayNorm_0025005<-SCnorm_runMAST(Data=cbind(mbayNorm_out_0025005$Bay_mat_list$`Group 1`,mbayNorm_out_0025005$Bay_mat_list$`Group 2`),NumCells = gr)



abayNorm_out_0025005<-bayNorm_sup(Data=RAW_DAT,BETA_vec = unlist(mbayNorm_out_0025005$BETA),PRIORS = mbayNorm_out_0025005$PRIORS_LIST,Conditions = CONDITION,S=10)

library(abind)
qq<-abind(abayNorm_out_0025005$Bay_array_list,along=2)
M_abayNorm_0025005<-SCnorm_runMAST3(Data=qq,NumCells=as.numeric(table(CONDITION)))


save(mbayNorm_out_0025005,M_mbayNorm_0025005,abayNorm_out_0025005,M_abayNorm_0025005,file="E:/RNAseqProject/SIMULATION/SIM_005_01/menabeta_0025005.RData")

###### 0.1 0.2###########
library(scran)
BETA_SCRAN_1_01<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1_01<-BETA_SCRAN_1_01/mean(BETA_SCRAN_1_01)*0.1
BETA_SCRAN_2_02<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2_02<-BETA_SCRAN_2_02/mean(BETA_SCRAN_2_02)*0.2



mbayNorm_out_0102<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1_01,BETA_SCRAN_2_02),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'LL')

gr<-table(CONDITION)
M_mbayNorm_0102<-SCnorm_runMAST(Data=cbind(mbayNorm_out_0102$Bay_mat_list$`Group 1`,mbayNorm_out_0102$Bay_mat_list$`Group 2`),NumCells = gr)
T_mbayNorm_0102<-wilcox_fun(cells=mbayNorm_out_0102$Bay_mat_list$`Group 1`,ctrls=mbayNorm_out_0102$Bay_mat_list$`Group 2`)


abayNorm_out_0102<-bayNorm_sup(Data=RAW_DAT,BETA_vec = unlist(mbayNorm_out_0102$BETA),PRIORS = mbayNorm_out_0102$PRIORS_LIST,Conditions = CONDITION,S=10)

library(abind)
qq<-abind(abayNorm_out_0102$Bay_array_list,along=2)
M_abayNorm_0102<-SCnorm_runMAST3(Data=qq,NumCells=as.numeric(table(CONDITION)))


save(mbayNorm_out_0102,M_mbayNorm_0102,T_mbayNorm_0102,abayNorm_out_0102,M_abayNorm_0102,file="E:/RNAseqProject/SIMULATION/SIM_005_01/menabeta_0102.RData")



######analysis#########
load("E:/RNAseqProject/SIMULATION/SIM_005_01/SIM_005_01.RData")
load("E:/RNAseqProject/SIMULATION/SIM_005_01/menabeta_0102.RData")
load("E:/RNAseqProject/SIMULATION/SIM_005_01/menabeta_0025005.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")




##ROC
MAST_RE_LIST<-list(bayNorm_10_00501=apply(M_bay10,1,median),bayNorm_10_0025005=apply(M_abayNorm_0025005,1,median),bayNorm_10_0102=apply(M_abayNorm_0102,1,median),bayNorm_mean_00501=M_mbayNorm$adjpval,bayNorm_mean_0025005=M_mbayNorm_0025005$adjpval,bayNorm_mean_0102=M_mbayNorm_0102$adjpval)



library(foreach)
TRUE_LABEL_reverse<-TRUE_LABEL
TRUE_LABEL_reverse[TRUE_LABEL_reverse==1]=3
TRUE_LABEL_reverse[TRUE_LABEL_reverse==0]=1
TRUE_LABEL_reverse[TRUE_LABEL_reverse==3]=0

Sele<-'MAST'
method_vec<-c('bayNorm_10_5%_10%:','bayNorm_10_2.5%_5%:','bayNorm_10_10%_20%:','bayNorm_mean_5%_10%:','bayNorm_mean_2.5%_5%:','bayNorm_mean_10%_20%:')
#col_vec<-seq(1,length(method_vec))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('bayNorm_10_5%_10%:','bayNorm_10_2.5%_5%:','bayNorm_10_10%_20%:','bayNorm_mean_5%_10%:','bayNorm_mean_2.5%_5%:','bayNorm_mean_10%_20%:','NULL1','NULL2')
col_vec<-cbbPalette[method_vec]


Input_re_list<-MAST_RE_LIST
mainn='Different normalization methods, MAST for DE detection'




auc_vec_M<-NULL
TRUE_LABEL_input<-TRUE_LABEL
table(TRUE_LABEL_input)

library(ROCR)
list_pref_M<-foreach(i=1:length(Input_re_list))%do%{
    pred_MAST <- prediction((Input_re_list[[i]]), TRUE_LABEL_reverse)
    perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
    
    auc_temp<-performance( pred_MAST, measure='auc' )
    auc_temp<-auc_temp@y.values[[1]]
    auc_vec_M<-c(auc_vec_M,auc_temp)
    return(perf_MAST)
}



#######Fig S24(b)#####

ROC_fun(list_pref=list_pref_M,vec_auc=auc_vec_M,method_vec=method_vec,col_vec=col_vec,MAIN=paste(''),cex=1,cex.axis=1,lwd=1.5,cex.lab=1,cex.legend=1,line=2)
abline(0,1,lty=2)
