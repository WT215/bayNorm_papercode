load("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/Torre_CV_default_divmeanbeta.RData")
load("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/Torre_GINI_default_divmeanbeta.RData")

load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_smFISH_meanBETA/CV_MEAN_data_default_divmeanbeta.RData")
load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_smFISH_meanBETA/GINI_data_default_divmeanbeta.RData")


library(ggplot2)
#paper
textsize<-10

#group meeting
textsize<-14
###Torre_CV#########

BAR_DAT_torrenew<-BAR_DAT
BAR_DAT_torrenew$CV<-log2(BAR_DAT_torrenew$CV/BAR_DAT$CV[BAR_DAT$Method=='smFISH'])


BAR_DAT_torrenew<-as.data.frame(BAR_DAT_torrenew)
BAR_DAT_torrenew$CV<-as.numeric(BAR_DAT_torrenew$CV)
BAR_DAT_torrenew$Method<-factor(BAR_DAT_torrenew$Method,levels=unique(BAR_DAT_torrenew$Method))
BAR_DAT_torrenew$Genename<-factor(BAR_DAT_torrenew$Genename,levels=unique(BAR_DAT_torrenew$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Torre_CV=ggplot(BAR_DAT_torrenew[BAR_DAT_torrenew$Method!='smFISH',],aes(x=Method,y=CV,col=Method,fill=Method))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of CV: scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Torre_CV



####Torre Gini######

BAR_out_Torre_GINI_v2<-BAR_out_Torre_GINI
BAR_out_Torre_GINI_v2$Gini<-log2(BAR_out_Torre_GINI_v2$Gini/BAR_out_Torre_GINI$Gini[BAR_out_Torre_GINI$Method=='smFISH'])


BAR_out_Torre_GINI_v2<-as.data.frame(BAR_out_Torre_GINI_v2)
BAR_out_Torre_GINI_v2$Gini<-as.numeric(BAR_out_Torre_GINI_v2$Gini)
BAR_out_Torre_GINI_v2$Method<-factor(BAR_out_Torre_GINI_v2$Method,levels=unique(BAR_out_Torre_GINI_v2$Method))
BAR_out_Torre_GINI_v2$Genename<-factor(BAR_out_Torre_GINI_v2$Genename,levels=unique(BAR_out_Torre_GINI_v2$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Torre_Gini=ggplot(BAR_out_Torre_GINI_v2[BAR_out_Torre_GINI_v2$Method!='smFISH',],aes(x=Method,y=Gini,col=Method,fill=Method))+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.65,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of Gini: scRNA-seq / smFISH")+
  
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Torre_Gini

####Torre mean##########


MeanExp_data_Torre_v2<-MeanExp_data_Torre
MeanExp_data_Torre_v2<-as.data.frame(MeanExp_data_Torre_v2)

MeanExp_data_Torre_v2$MeanExp<-as.numeric(as.character(MeanExp_data_Torre_v2$MeanExp))
MeanExp_data_Torre_v2$Gene<-factor(MeanExp_data_Torre_v2$Gene,levels=unique(MeanExp_data_Torre_v2$Gene))
MeanExp_data_Torre_v2$NormMethods<-factor(MeanExp_data_Torre_v2$NormMethods,levels=c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA'))

smFISH_mean<-MeanExp_data_Torre_v2$MeanExp[MeanExp_data_Torre_v2$NormMethods=='smFISH']
MeanExp_data_Torre_v2$MeanExp<-log2(MeanExp_data_Torre_v2$MeanExp/smFISH_mean)

Torre_mean=ggplot(MeanExp_data_Torre_v2[MeanExp_data_Torre_v2$NormMethods!='smFISH',],aes(x=NormMethods,y=MeanExp,col=NormMethods,fill=NormMethods))+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of mean: scRNA-seq / smFISH")+
  
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Torre_mean


#########Grun serum#########

###Grun S_CV#########

BAR_DAT_S_v2<-BAR_DAT_S
BAR_DAT_S_v2$CV<-log2(BAR_DAT_S_v2$CV/BAR_DAT_S$CV[BAR_DAT_S$Method=='smFISH'])


BAR_DAT_S_v2<-as.data.frame(BAR_DAT_S_v2)
BAR_DAT_S_v2$CV<-as.numeric(BAR_DAT_S_v2$CV)
BAR_DAT_S_v2$Method<-factor(BAR_DAT_S_v2$Method,levels=unique(BAR_DAT_S_v2$Method))
BAR_DAT_S_v2$Genename<-factor(BAR_DAT_S_v2$Genename,levels=unique(BAR_DAT_S_v2$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

GrunS_CV=ggplot(BAR_DAT_S_v2[BAR_DAT_S_v2$Method!='smFISH',],aes(x=Method,y=CV,col=Method,fill=Method))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
GrunS_CV



####Grun S Gini######

BAR_DAT_S_GINI_v2<-BAR_DAT_S_GINI
BAR_DAT_S_GINI_v2$Gini<-log2(BAR_DAT_S_GINI_v2$Gini/BAR_DAT_S_GINI$Gini[BAR_DAT_S_GINI$Method=='smFISH'])


BAR_DAT_S_GINI_v2<-as.data.frame(BAR_DAT_S_GINI_v2)
BAR_DAT_S_GINI_v2$Gini<-as.numeric(BAR_DAT_S_GINI_v2$Gini)
BAR_DAT_S_GINI_v2$Method<-factor(BAR_DAT_S_GINI_v2$Method,levels=unique(BAR_DAT_S_GINI_v2$Method))
BAR_DAT_S_GINI_v2$Genename<-factor(BAR_DAT_S_GINI_v2$Genename,levels=unique(BAR_DAT_S_GINI_v2$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

GrunS_Gini=ggplot(BAR_DAT_S_GINI_v2[BAR_DAT_S_GINI_v2$Method!='smFISH',],aes(x=Method,y=Gini,col=Method,fill=Method))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
GrunS_Gini

####Grun S mean##########


MeanExp_data_S_v2<-MeanExp_data_serum
MeanExp_data_S_v2<-as.data.frame(MeanExp_data_S_v2)

MeanExp_data_S_v2$MeanExp<-as.numeric(as.character(MeanExp_data_S_v2$MeanExp))
MeanExp_data_S_v2$Gene<-factor(MeanExp_data_S_v2$Gene,levels=unique(MeanExp_data_S_v2$Gene))
MeanExp_data_S_v2$NormMethods<-factor(MeanExp_data_S_v2$NormMethods,levels=c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA'))

smFISH_mean<-MeanExp_data_S_v2$MeanExp[MeanExp_data_S_v2$NormMethods=='smFISH']
MeanExp_data_S_v2$MeanExp<-log2(MeanExp_data_S_v2$MeanExp/smFISH_mean)

GrunS_mean=ggplot(MeanExp_data_S_v2[MeanExp_data_S_v2$NormMethods!='smFISH',],aes(x=NormMethods,y=MeanExp,col=NormMethods,fill=NormMethods))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
GrunS_mean



#########Grun 2i#########
###Grun 2i_CV#########

BAR_DAT_2i_v2<-BAR_DAT_2i
BAR_DAT_2i_v2$CV<-log2(BAR_DAT_2i_v2$CV/BAR_DAT_2i$CV[BAR_DAT_2i$Method=='smFISH'])


BAR_DAT_2i_v2<-as.data.frame(BAR_DAT_2i_v2)
BAR_DAT_2i_v2$CV<-as.numeric(BAR_DAT_2i_v2$CV)
BAR_DAT_2i_v2$Method<-factor(BAR_DAT_2i_v2$Method,levels=unique(BAR_DAT_2i_v2$Method))
BAR_DAT_2i_v2$Genename<-factor(BAR_DAT_2i_v2$Genename,levels=unique(BAR_DAT_2i_v2$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#combine 2i and serum results for fig 2####
inputdata<-rbind(BAR_DAT_2i_v2,BAR_DAT_S_v2)

Grun2i_CV=ggplot(inputdata[inputdata$Method!='smFISH',],aes(x=Method,y=CV,col=Method,fill=Method))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of CV: scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Grun2i_CV



####Grun 2i Gini######

BAR_DAT_2i_GINI_v2<-BAR_DAT_2i_GINI
BAR_DAT_2i_GINI_v2$Gini<-log2(BAR_DAT_2i_GINI_v2$Gini/BAR_DAT_2i_GINI$Gini[BAR_DAT_2i_GINI$Method=='smFISH'])


BAR_DAT_2i_GINI_v2<-as.data.frame(BAR_DAT_2i_GINI_v2)
BAR_DAT_2i_GINI_v2$Gini<-as.numeric(BAR_DAT_2i_GINI_v2$Gini)
BAR_DAT_2i_GINI_v2$Method<-factor(BAR_DAT_2i_GINI_v2$Method,levels=unique(BAR_DAT_2i_GINI_v2$Method))
BAR_DAT_2i_GINI_v2$Genename<-factor(BAR_DAT_2i_GINI_v2$Genename,levels=unique(BAR_DAT_2i_GINI_v2$Genename))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


inputdata<-rbind(BAR_DAT_2i_GINI_v2,BAR_DAT_S_GINI_v2)

Grun2i_Gini=ggplot(inputdata[inputdata$Method!='smFISH',],aes(x=Method,y=Gini,col=Method,fill=Method))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of Gini: scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Grun2i_Gini

####Grun 2i mean##########


MeanExp_data_2i_v2<-MeanExp_data_2i
MeanExp_data_2i_v2<-as.data.frame(MeanExp_data_2i_v2)

MeanExp_data_2i_v2$MeanExp<-as.numeric(as.character(MeanExp_data_2i_v2$MeanExp))
MeanExp_data_2i_v2$Gene<-factor(MeanExp_data_2i_v2$Gene,levels=unique(MeanExp_data_2i_v2$Gene))
MeanExp_data_2i_v2$NormMethods<-factor(MeanExp_data_2i_v2$NormMethods,levels=c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA'))

smFISH_mean<-MeanExp_data_2i_v2$MeanExp[MeanExp_data_2i_v2$NormMethods=='smFISH']
MeanExp_data_2i_v2$MeanExp<-log2(MeanExp_data_2i_v2$MeanExp/smFISH_mean)

inputdata<-rbind(MeanExp_data_2i_v2,MeanExp_data_S_v2)

Grun2i_mean=ggplot(inputdata[inputdata$NormMethods!='smFISH',],aes(x=NormMethods,y=MeanExp,col=NormMethods,fill=NormMethods))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
  scale_fill_manual(values=cbbPalette[-1])+
  labs(x='',y="log2 ratio of mean: scRNA-seq / smFISH")+
  geom_boxplot(position="identity",col="black",fill=NA,width=0.2,outlier.colour =NA)+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-1,1),lty=2)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
Grun2i_mean


qq<-ggarrange(Grun2i_mean + theme(legend.position="none"),Grun2i_CV + theme(legend.position="none"),Grun2i_Gini+ theme(legend.position="none"),ncol=3,nrow=1,common.legend = TRUE, legend="none")
qq
ggsave(filename ="C:/Users/Wenhao/Desktop/Grun_default.pdf",plot = qq,width=16,height=5.5)


######FIG 2 #########

library(gridExtra)
library(ggpubr)
library(cowplot)

# qq<-plot_grid(Grun2i_mean + theme(legend.position="none"),Grun2i_CV + theme(legend.position="none"),Grun2i_Gini+ theme(legend.position="none"), Torre_mean + theme(legend.position="none"),Torre_CV + theme(legend.position="none"),Torre_Gini+ theme(legend.position="none"),ncol=3,nrow=2)+ draw_grob(get_legend(Grun2i_mean), (1.5)/3, 0, 1/3, 0.5)
# qq


qq<-ggarrange(Grun2i_mean + theme(legend.position="none"),Grun2i_CV + theme(legend.position="none"),Grun2i_Gini+ theme(legend.position="none"), Torre_mean + theme(legend.position="none"),Torre_CV + theme(legend.position="none"),Torre_Gini+ theme(legend.position="none"),ncol=3,nrow=2,common.legend = TRUE, legend="none")
qq
ggsave(filename ="E:/RNAseqProject/Illustrator_bayNorm/newFIG/Grun_Torre_default_divmeanbeta.pdf",plot = qq,width=16,height=11)
#ggsave(filename ="E:/IC-PhD/GroupMeeting/22062018/Figure_22062018/Grun_Torre_default.pdf",plot = qq,width=19,height=13.5)
