#Klein based simulation
load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005/SIM_noDE_01_005.RData")
load("E:/RNAseqProject/SIMULATION/DE_GG_explore/GG_SIM_noDE_01_005.RData")




length(which(M_LL_mbayNorm_out$adjpval<0.05))

FPR<-c(length(which(apply(M_LL_bayNorm_out,1,median)<0.05))/9999,length(which(M_LL_mbayNorm_out$adjpval<0.05))/9999,length(which(apply(M_bayNorm_a10,1,median)<0.05))/9999,length(which(M_mbayNorm$adjpval<0.05))/9999)
bardata<-data.frame(FPR,c('Local','Local','Global','Global'),c('3D array (10 samples)','mean of posterior','3D array (10 samples)','mean of posterior'))
colnames(bardata)<-c('FPR','type of procedures','type of output')




library(ggplot2)
textsize<-14
FPR_bar<-ggplot(data=bardata, aes(x=bardata$`type of output`, y=bardata$FPR, fill=bardata$`type of procedures`)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=round(FPR,4)), vjust=0, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "Type of output from bayNorm",y='False positive rates',fill='Prior parameters')+ggtitle("") +
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
FPR_bar

# g <- ggplot_build(FPR_bar)
# unique(g$data[[1]]["fill"])

####FC

DATA_list<-list(Global=cbind(mbayNorm_out$Bay_mat_list$`Group 1`,mbayNorm_out$Bay_mat_list$`Group 2`),Local=cbind(LL_mbayNorm_out$Bay_mat_list$`Group 1`,LL_mbayNorm_out$Bay_mat_list$`Group 2`))

source("E:/RNAseqProject/Bacher__SCnorm_2016/FC_fun.R")
FC_out<-FC_fun(Inputdat=RAW_DAT,CONDITION=CONDITION,DATA_list,textsize=14,legend.key.size=1,colourval = c('#F8766D','#00BFC4'))
FC_out



##Begin plot for Fig S20 (a)-(b)####

library(gridExtra)
library(ggpubr)
library(cowplot)

qq<-plot_grid(FPR_bar + theme(legend.position=c(0.01,0.85)),FC_out + theme(legend.position="none"),ncol=2,nrow=1)
# qq<-plot_grid(FPR_bar + theme(legend.position=c(0.01,0.85)),FC_out + theme(legend.position="none"),AUC_bar_a + theme(legend.position="none"),AUC_bar_m+theme(legend.position="none"),ncol=2,nrow=2)
# qq
