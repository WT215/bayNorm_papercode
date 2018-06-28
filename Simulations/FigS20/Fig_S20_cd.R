
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



load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_1.RData")
(auc_aGG_0101<-auc_GG_bayNorm_out)
(auc_mGG_0101<-auc_GG_mbayNorm_out)

load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_01_005.RData")
(auc_aGG_01005<-auc_GG_bayNorm_out)
(auc_mGG_01005<-auc_GG_mbayNorm_out)


load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_005_01.RData")
(auc_aGG_00501<-auc_GG_bayNorm_out)
(auc_mGG_00501<-auc_GG_mbayNorm_out)

load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_005_005.RData")
(auc_aGG_005005<-auc_GG_bayNorm_out)
(auc_mGG_005005<-auc_GG_mbayNorm_out)



save(auc_aGG_0101,auc_mGG_0101,auc_aGG_01005,auc_mGG_01005,auc_aGG_00501,auc_mGG_00501,auc_aGG_005005,auc_mGG_005005,file="E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_aucs.RData")

names(S1_M_list2)
###LL#######
(auc_aLL_0101<-auc_fun(pval=S1_M_list2$bayNorm_10,TRUELABEL = TRUE_LABEL_reverse_S1))
(auc_mLL_0101<-auc_fun(pval=S1_M_list2$bayNorm_mean,TRUELABEL = TRUE_LABEL_reverse_S1))

(auc_aLL_01005<-auc_fun(pval=S4_M_list2$bayNorm_10,TRUELABEL = TRUE_LABEL_reverse_S4))
(auc_mLL_01005<-auc_fun(pval=S4_M_list2$bayNorm_mean,TRUELABEL = TRUE_LABEL_reverse_S4))


(auc_aLL_00501<-auc_fun(pval=S2_M_list2$bayNorm_10,TRUELABEL = TRUE_LABEL_reverse_S2))
(auc_mLL_00501<-auc_fun(pval=S2_M_list2$bayNorm_mean,TRUELABEL = TRUE_LABEL_reverse_S2))

(auc_aLL_005005<-auc_fun(pval=S3_M_list2$bayNorm_10,TRUELABEL = TRUE_LABEL_reverse_S3))
(auc_mLL_005005<-auc_fun(pval=S3_M_list2$bayNorm_mean,TRUELABEL = TRUE_LABEL_reverse_S3))

save(auc_aLL_0101,auc_mLL_0101,auc_aLL_01005,auc_mLL_01005,auc_aLL_00501,auc_mLL_00501,auc_aLL_005005,auc_mLL_005005,file="E:/RNAseqProject/SIMULATION/DE_GG_explore/LL_aucs.RData")



####begin analysis##########
load("E:/RNAseqProject/SIMULATION/DE_GG_explore/LL_aucs.RData")
load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_aucs.RData")

varnames<-ls()

library(foreach)
library(stringr)
str_sub(varnames,9,-1)

aucdata<-foreach(i=1:length(varnames),.combine=rbind)%do%{
    ee<-c(varnames[i],get(varnames[i]),str_sub(varnames[i],6,7),str_sub(varnames[i],5,5),str_sub(varnames[i],9,-1))
    return(ee)
}
colnames(aucdata)<-c('varnames','AUC','type of procedure','type of output','SIM cases')
class(aucdata)
aucdata<-as.data.frame(aucdata)

aucdata$AUC<-as.numeric(as.character(aucdata$AUC))
aucdata$`SIM cases`<-rep(c('SIM DE IV','SIM DE II','SIM DE III','SIM DE I'),4)
aucdata$`type of procedure`<-rep(c(rep('Global',4),rep('Local',4)),2)
aucdata$`type of output`

library(ggplot2)
textsize<-14
AUC_bar_m<-ggplot(data=aucdata[which(aucdata$`type of output`=='m'),], aes(x=`SIM cases`, y=AUC, fill=`type of procedure`)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=round(AUC,4)), vjust=0, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "SIM cases",y='AUC',fill='Prior parameters')+ggtitle("") +
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
AUC_bar_m



library(ggplot2)
textsize<-14
AUC_bar_a<-ggplot(data=aucdata[which(aucdata$`type of output`=='a'),], aes(x=`SIM cases`, y=AUC, fill=`type of procedure`)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=round(AUC,4)), vjust=0, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "SIM cases",y='AUC',fill='Prior parameters')+ggtitle("") +
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
AUC_bar_a



##Fig S20 (c)-(d)
source(file="E:/RNAseqProject/MANY_SAVE_PATH.R")
library(gridExtra)
library(ggpubr)
library(cowplot)
qq<-plot_grid(AUC_bar_a,AUC_bar_m + theme(legend.position="none"),ncol=1,nrow=2)
qq

