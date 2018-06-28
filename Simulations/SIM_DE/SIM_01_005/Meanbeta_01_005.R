load("E:/RNAseqProject/SIMULATION/SIM_01_005/SIM_01_005.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")
library(bayNorm)

library(scran)
BETA_SCRAN_1_005<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1_005<-BETA_SCRAN_1_005/mean(BETA_SCRAN_1_005)*0.05
BETA_SCRAN_2_0025<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2_0025<-BETA_SCRAN_2_0025/mean(BETA_SCRAN_2_0025)*0.025



#####0.05 0.025#########

mbayNorm_out_0050025<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1_005,BETA_SCRAN_2_0025),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'LL')

gr<-table(CONDITION)
M_mbayNorm_0050025<-SCnorm_runMAST(Data=cbind(mbayNorm_out_0050025$Bay_mat_list$`Group 1`,mbayNorm_out_0050025$Bay_mat_list$`Group 2`),NumCells = gr)


abayNorm_out_0050025<-bayNorm_sup(Data=RAW_DAT,BETA_vec = unlist(mbayNorm_out_0050025$BETA),PRIORS = mbayNorm_out_0050025$PRIORS_LIST,Conditions = CONDITION,S=10)

library(abind)
qq<-abind(abayNorm_out_0050025$Bay_array_list,along=2)
M_abayNorm_0050025<-SCnorm_runMAST3(Data=qq,NumCells=as.numeric(table(CONDITION)))


save(mbayNorm_out_0050025,M_mbayNorm_0050025,T_mbayNorm_0050025,abayNorm_out_0050025,M_abayNorm_0050025,file="E:/RNAseqProject/SIMULATION/SIM_01_005/menabeta_0050025.RData")

######0.2 0.1###########
library(scran)
BETA_SCRAN_1_02<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1_02<-BETA_SCRAN_1_02/mean(BETA_SCRAN_1_02)*0.2
BETA_SCRAN_2_01<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2_01<-BETA_SCRAN_2_01/mean(BETA_SCRAN_2_01)*0.1



mbayNorm_out_0201<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1_02,BETA_SCRAN_2_01),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'LL')

gr<-table(CONDITION)
M_mbayNorm_0201<-SCnorm_runMAST(Data=cbind(mbayNorm_out_0201$Bay_mat_list$`Group 1`,mbayNorm_out_0201$Bay_mat_list$`Group 2`),NumCells = gr)


abayNorm_out_0201<-bayNorm_sup(Data=RAW_DAT,BETA_vec = unlist(mbayNorm_out_0201$BETA),PRIORS = mbayNorm_out_0201$PRIORS_LIST,Conditions = CONDITION,S=10)

library(abind)
qq<-abind(abayNorm_out_0201$Bay_array_list,along=2)
M_abayNorm_0201<-SCnorm_runMAST3(Data=qq,NumCells=as.numeric(table(CONDITION)))


save(mbayNorm_out_0201,M_mbayNorm_0201,T_mbayNorm_0201,abayNorm_out_0201,M_abayNorm_0201,file="E:/RNAseqProject/SIMULATION/SIM_01_005/menabeta_0201.RData")
load("E:/RNAseqProject/SIMULATION/SIM_01_005/Meanbetaexplore/menabeta_0201.RData")


######analysis#########
load("E:/RNAseqProject/SIMULATION/SIM_01_005/SIM_01_005.RData")
load("E:/RNAseqProject/SIMULATION/SIM_01_005/menabeta_0201.RData")
load("E:/RNAseqProject/SIMULATION/SIM_01_005/menabeta_0050025.RData")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")




##ROC
MAST_RE_LIST<-list(bayNorm_10_01005=apply(M_bay10,1,median),bayNorm_10_0050025=apply(M_abayNorm_0050025,1,median),bayNorm_10_0201=apply(M_abayNorm_0201,1,median),bayNorm_mean_01005=M_mbayNorm$adjpval,bayNorm_mean_0050025=M_mbayNorm_0050025$adjpval,bayNorm_mean_0201=M_mbayNorm_0201$adjpval)



library(foreach)
TRUE_LABEL_reverse<-TRUE_LABEL
TRUE_LABEL_reverse[TRUE_LABEL_reverse==1]=3
TRUE_LABEL_reverse[TRUE_LABEL_reverse==0]=1
TRUE_LABEL_reverse[TRUE_LABEL_reverse==3]=0

Sele<-'MAST'
method_vec<-c('bayNorm_10_10%_5%:','bayNorm_10_5%_2.5%:','bayNorm_10_20%_10%:','bayNorm_mean_10%_5%:','bayNorm_mean_5%_2.5%:','bayNorm_mean_20%_10%:')
#col_vec<-seq(1,length(method_vec))


#col_vec<-seq(1,length(method_vec))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('bayNorm_10_10%_5%:','bayNorm_10_5%_2.5%:','bayNorm_10_20%_10%:','bayNorm_mean_10%_5%:','bayNorm_mean_5%_2.5%:','bayNorm_mean_20%_10%:','NULL1','NULL2')
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




#######Fig S24 (c)#####
ROC_fun(list_pref=list_pref_M,vec_auc=auc_vec_M,method_vec=method_vec,col_vec=col_vec,MAIN='',cex=1,cex.axis=1,lwd=1.5,cex.lab=1,cex.legend=1,line=2)
abline(0,1,lty=2)

