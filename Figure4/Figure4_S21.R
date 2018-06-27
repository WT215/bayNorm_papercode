load(file='Tung_norms.RData')



MAST_RE_LIST<-list(MAST_bay=apply(M_bay_mat,1,median),MAST_scnorm=M_scnorm$adjpval,MAST_scimpute=MAST_scimpute$adjpval,MAST_RB=MAST_RB$adjpval,SAVER=apply(MAST_saver_array_mat,1,median),MAGIC=M_MAGIC_tr$adjpval,DCA=M_DCA$adjpval)
Sele<-'MAST'
method_vec<-c('bayNorm','SCnorm','scImpute','Scaling','SAVER','MAGIC','DCA')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
col_vec<-cbbPalette[method_vec]
Input_re_list<-MAST_RE_LIST
mainn='DE detection method: MAST'

TRUE_LABEL_input<-DE_TRUE_LABEL

TRUE_LABEL_input[TRUE_LABEL_input==0]=3
TRUE_LABEL_input[TRUE_LABEL_input==1]=0
TRUE_LABEL_input[TRUE_LABEL_input==3]=1


table(DE_TRUE_LABEL)
length(DE_TRUE_LABEL)

auc_vec<-NULL
library(ROCR)
list_pref<-foreach(i=1:length(Input_re_list))%do%{
    pred_MAST <- prediction(Input_re_list[[i]], TRUE_LABEL_input)
    perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
    
    auc_temp<-performance( pred_MAST, measure='auc' )
    auc_temp<-auc_temp@y.values[[1]]
    auc_vec<-c(auc_vec,auc_temp)
    return(perf_MAST)
}


#AUC value in Figure 4(c)########
ROC_fun(list_pref=list_pref,vec_auc=auc_vec,method_vec=method_vec,col_vec=col_vec,MAIN=mainn,cex=1,cex.axis=1,lwd=1.5,cex.lab=1,cex.legend=1,line=2)
abline(0,1,lty=2)

BATCH_ALL<-c(rep(1,dim(N1_1_DAT)[2]),rep(3,dim(N1_3_DAT)[2]),rep(1,dim(N2_1_DAT)[2]),rep(2,dim(N2_2_DAT)[2]),rep(3,dim(N2_3_DAT)[2]),rep(1,dim(N3_1_DAT)[2]),rep(2,dim(N3_2_DAT)[2]),rep(3,dim(N3_3_DAT)[2]))
names(BATCH_ALL)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))


#Figure 4(b): PCA: bayNorm####
library(abind)
#PCA: mean version's bayNorm
#pca_try <- prcomp(t(cbind(bayNorm_mean_N1$Bay_mat,bayNorm_mean_N2$Bay_mat,bayNorm_mean_N3$Bay_mat)),scale=T)
source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")
textsize=6
pointsize=0.5
legendpointsize=1.5
legend_key_size<-1

#PCA: one sample from bayNorm
pca_array <- prcomp(t(cbind(Bay_1[,,1],Bay_2[,,1],Bay_3[,,1])),scale=T)
BATCH_array<-PCA_FUN(pca_try=pca_array,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,legend_key_size=legend_key_size,TITLE='bayNorm (1 sample)')
BATCH_array

#Figure 4(a): PCA: Scaling method####
#Scaling method
pca_try_RB <- prcomp(t(cbind(RB_N1,RB_N2,RB_N3)),scale=T)
BATCH_RB<-PCA_FUN(pca_try=pca_try_RB,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='Scaling')


#Fig S25######
pca_try_SAVER <- prcomp(t(cbind(saver_N1$estimate,saver_N2$estimate,saver_N3$estimate)),scale=T)
pca_try_SCnorm <- prcomp(t(scnorm_norm),scale=T)
pca_try_scImpute <- prcomp(t(cbind(scImpute_N1_DAT,scImpute_N2_DAT,scImpute_N3_DAT)),scale=T)
pca_try_MAGIC <- prcomp(t(MAGIC_TUNG),scale=T)
pca_try_DCA <- prcomp(t(DCA_tungall),scale=T)




source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")
textsize=14
pointsize=1
legendpointsize=1
BATCH_SAVER<-PCA_FUN(pca_try=pca_try_SAVER,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='SAVER')
BATCH_SCnorm<-PCA_FUN(pca_try=pca_try_SCnorm,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='SCnorm')
BATCH_scImpute<-PCA_FUN(pca_try=pca_try_scImpute,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='scImpute')
BATCH_MAGIC<-PCA_FUN(pca_try=pca_try_MAGIC,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='MAGIC')

BATCH_DCA<-PCA_FUN(pca_try=pca_try_DCA,LABEL_INDIVIDUAL =LABEL_INDIVIDUAL ,LABEL_REP = LABEL_REP,textsize=textsize,pointsize = pointsize,legendpointsize=legendpointsize,TITLE='MAGIC')

line<-2.2
line_title<-0.5
cex.lab<-1
cex<-1
par(mfrow=c(3,2))
PCA_FUN2(pca_try=pca_try_SAVER ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='SAVER',legend.x='topleft')
PCA_FUN2(pca_try=pca_try_SCnorm ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='SCnorm',legend.x='topleft')
PCA_FUN2(pca_try=pca_try_scImpute ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='scImpute',legend.x='topleft')
PCA_FUN2(pca_try=pca_try_MAGIC ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='MAGIC',legend.x='topleft')
PCA_FUN2(pca_try=pca_try_DCA ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='DCA',legend.x='topleft')






