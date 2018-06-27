load('Load_Islam.RData')
load('Islam_many_normalizations.RData')



MAST_RE_LIST<-list(MAST_bay=apply(MAST_bay_mat,1,median),MAST_SCnorm=MAST_RE_LIST$MAST_SCnorm,MAST_scImpute=M_scImpute_Islam$adjpval,Scaling=MAST_RB$adjpval,MAGIC=M_magic_islam$adjpval,DCA=M_DCA$adjpval)



library(foreach)
source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

mmme<-'MAST'


Input_re_list<-MAST_RE_LIST





TRUE_LABEL_input<-TRUE_LABEL


TRUE_LABEL_input[TRUE_LABEL_input==0]=3
TRUE_LABEL_input[TRUE_LABEL_input==1]=0
TRUE_LABEL_input[TRUE_LABEL_input==3]=1


library(ROCR)
list_pref<-foreach(i=1:length(Input_re_list))%do%{
    pred_MAST <- prediction(Input_re_list[[i]], TRUE_LABEL_input)
    perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
    
    auc_temp<-performance( pred_MAST, measure='auc' )
    auc_temp<-auc_temp@y.values[[1]]
    auc_vec<-c(auc_vec,auc_temp)
    return(perf_MAST)
}

method_vec<-c('bayNorm','SCnorm','scImpute','Scaling','MAGIC','DCA')


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('bayNorm_10','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
col_vec<-cbbPalette[method_vec]

source("E:/RNAseqProject/MANY_DE_FUN.R")

ROC_fun(list_pref=list_pref,vec_auc=auc_vec,method_vec=method_vec,col_vec=col_vec,MAIN=paste('DE detection method:', mmme),cex=10,lwd=3,cex.lab=10,cex.axis=10,cex.legend=1.5)
abline(0,1)
