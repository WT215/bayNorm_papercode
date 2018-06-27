load('Grun_serum_norms.RData')
load('Grun_serum_2i.RData')
load("Grun_2014_RAW.RData")
load('Grun_2014_RAW_serum.RData')
load("smFISH_norm_load.RData")

source("E:/RNAseqProject/Density_fun.r")

SCnorm_2i<-SCnorm_2i_ori@metadata$NormalizedData
SCnorm_serum<-SCnorm_serum_ori@metadata$NormalizedData



bay_2i_input<-bayNorm_SC_2i$Bay_array
bay_serum_input<-bayNorm_SC_serum$Bay_array


#Scale normalized data by mean BETA used in bayNorm, for fair comparison
MAGIC_2i<-MAGIC_2i/mean(bayNorm_SC_2i$BETA)
DCA_2i<-DCA_2i/mean(bayNorm_SC_2i$BETA)
scImpute_SC_2i<-scImpute_SC_2i/mean(bayNorm_SC_2i$BETA)
SAVER_SC_2i_array<-SAVER_SC_2i_array/mean(bayNorm_SC_2i$BETA)
SCnorm_2i<-SCnorm_2i/mean(bayNorm_SC_2i$BETA)
MAGIC_2i<-MAGIC_2i/mean(bayNorm_SC_2i$BETA)

DCA_serum<-DCA_serum/mean(bayNorm_SC_serum$BETA)
scImpute_SC_serum<-scImpute_SC_serum/mean(bayNorm_SC_serum$BETA)
SAVER_SC_serum_array<-SAVER_SC_serum_array/mean(bayNorm_SC_serum$BETA)
SCnorm_serum<-SCnorm_serum/mean(bayNorm_SC_serum$BETA)
MAGIC_serum<-MAGIC_serum/mean(bayNorm_SC_serum$BETA)

#Prepare smFISH data into list for convenience (2i)
LIST_2i<-Grun_2i_norm_list
Pou5f1_2i<-c(LIST_2i$Pou5f1_Klf4_2i[,1])
Fish_2iplot_list<-list(LIST_2i$Sohlh2_Hormad1_2i[,1],LIST_2i$Notch1_Pou5f1_2i[,1],LIST_2i$Stag3_Gli2_2i[,2],LIST_2i$Stag3_Gli2_2i[,1],LIST_2i$Tpx2_Pcna_J12i[,1],LIST_2i$Pou5f1_Klf4_2i[,2],LIST_2i$Tpx2_Pcna_J12i[,2],LIST_2i$Sox2_Pou5f1_2i[,1],Pou5f1_2i)
names(Fish_2iplot_list)<-Fig3Gene


#Prepare smFISH data into list for convenience (serum)
LIST_serum<-Grun_S_norm_list
Pou5f1_serum<-c(LIST_serum$Pou5f1_Klf4_J1S[,1])
Fish_serumplot_list<-list(LIST_serum$Sohlh2_Hormad1_S[,1],LIST_serum$Pou5f1_Notch1_S[,2],LIST_serum$Stag3_Gli2_S[,2],LIST_serum$Stag3_Gli2_S[,1],LIST_serum$Tpx2_Pcna_J1S[,1],LIST_serum$Pou5f1_Klf4_J1S[,2],LIST_serum$Tpx2_Pcna_J1S[,2],LIST_serum$Sox2_j1s[,1],Pou5f1_serum)
names(Fish_serumplot_list)<-Fig3Gene



#Preparing for the density plots:
########begin double for loops######

CNF<-function(x,numb=500){
  #qq<-x/mean(x)*numb
  #qq<-x/mean(x)
  qq<-x
  return(qq)
}

bw=10
gg_list_2i<-list()
gg_list_S<-list()
DenDat_list_2i<-list()
DenDat_list_S<-list()
ks_dist_2i<-list()
ks_dist_S<-list()



for(qq1 in seq(1,2)){
  
  typp<-c('2i','S')[qq1]
  print(typp)
  for(qq2 in seq(1,9)){
    
    Gg<-Fig3Gene[qq2]
    print(Gg)
    
    
    if(typp=="2i"){
      cccho<-which(names(Fish_2iplot_list)==Gg)
      fish_input<-CNF(Fish_2iplot_list[[cccho]],numb=nummm)
      nummm<-500
      bay_input<-CNF(as.vector(bay_2i_input[Gg,,]),numb=nummm)
      #SCnorm_input<-CNF(SCnorm_2i[Gg,],numb=nummm)
      scran_input<-CNF(scran_2i[Gg,],numb=nummm)
      #RPM_input<-CNF(RPM_2i[Gg,],numb=nummm)
      # TMM_input<-CNF(TMM_2i[Gg,],numb=nummm)
      #DESeq_input<-CNF(DESeq_2i[Gg,],numb=nummm)
      #SAVER_input<-CNF(SAVER_SC_2i$estimate[Gg,],numb=nummm)
      SAVER_input<-CNF(as.vector(SAVER_SC_2i_array[Gg,,]),numb=nummm)
      Raw_BETA_2i<-t(t(CountsUMI_SC_2i)/bayNorm_SC_2i$BETA)
      Raw_BETA_input<-CNF(Raw_BETA_2i[Gg,],numb=nummm)
      Raw_input<-CNF(CountsUMI_SC_2i[Gg,],numb=nummm)
      scImpute_input<-CNF(scImpute_SC_2i[Gg,],numb=nummm)
      
      
      ks_dist_2i[[qq2]]<-c(
        ks.test(bay_input,fish_input)$statistic,
        ks.test(Raw_BETA_input,fish_input)$statistic,
        ks.test(SAVER_input,fish_input)$statistic,
        ks.test(scImpute_input,fish_input)$statistic)
      
      
      inputlist<-list(smFISH=fish_input,bayNorm=bay_input,SAVER=SAVER_input,Scaling=Raw_BETA_input,Raw=Raw_input)
      #inputlist<-list(smFISH=fish_input,bayNorm=bay_input,Raw=Raw_input,Scaling=Raw_BETA_input)
      
      
    }else if(typp=="S"){
      cccho<-which(names(Fish_serumplot_list)==Gg)
      fish_input<-CNF(Fish_serumplot_list[[cccho]],numb=nummm)
      nummm<-500
      bay_input<-CNF(as.vector(bay_serum_input[Gg,,]),numb=nummm)
      SCnorm_input<-CNF(SCnorm_serum[Gg,],numb=nummm)
      scran_input<-CNF(scran_serum[Gg,],numb=nummm)
      # RPM_input<-CNF(RPM_serum[Gg,],numb=nummm)
      # TMM_input<-CNF(TMM_serum[Gg,],numb=nummm)
      # DESeq_input<-CNF(DESeq_serum[Gg,],numb=nummm)
      #SAVER_input<-CNF(SAVER_SC_serum$estimate[Gg,],numb=nummm)
      SAVER_input<-CNF(as.vector(SAVER_SC_serum_array[Gg,,]),numb=nummm)
      Raw_BETA_serum<-t(t(CountsUMI_SC_serum)/bayNorm_SC_serum$BETA)
      Raw_BETA_input<-CNF(Raw_BETA_serum[Gg,],numb=nummm)
      Raw_input<-CNF(CountsUMI_SC_serum[Gg,],numb=nummm)
      scImpute_input<-CNF(scImpute_SC_serum[Gg,],numb=nummm)
      
      
      ks_dist_S[[qq2]]<-c(
        ks.test(bay_input,fish_input)$statistic,
        ks.test(Raw_BETA_input,fish_input)$statistic,
        ks.test(SAVER_input,fish_input)$statistic,
        ks.test(scImpute_input,fish_input)$statistic)
      
      
      inputlist<-list(smFISH=fish_input,bayNorm=bay_input,SAVER=SAVER_input,Scaling=Raw_BETA_input,Raw=Raw_input)
      #inputlist<-list(smFISH=fish_input,bayNorm=bay_input,Raw=Raw_input,Scaling=Raw_BETA_input)
      
      
    }
    
    ########begin creating data for plotting density####
    
    #methodd<-c('smFISH','bayNorm','Scaling','SAVER','scImpute')
    #methodd<-c('smFISH','bayNorm','Scaling','Raw')
    methodd<-names(inputlist)
    
    
    library(foreach)
    DenDat<-foreach(i=1:length(inputlist),.combine=rbind)%do%{
      aa<-inputlist[[i]]
      bb<-cbind(aa,rep(methodd[i],length(aa)),rep(i,length(aa)),rep(Gg,length(aa)),rep(typp,length(aa)))
      return(bb)
    }
    colnames(DenDat)<-c('Normalized count','Normalization method','Colour','Genename','Condition')
    
    
    DenDat<-as.data.frame(DenDat)
    DenDat$`Normalized count`<-as.numeric(as.character(DenDat$`Normalized count`))
    DenDat$`Normalization method`<-factor(DenDat$`Normalization method`,levels=unique(DenDat$`Normalization method`))
    
    if(typp=='2i'){
      DenDat_list_2i[[qq2]]<-DenDat
      
    }else if(typp=='S'){
      DenDat_list_S[[qq2]]<-DenDat
    }
    
    
    
    filenn<-paste(typp,'_',Gg,sep='')
    library(ggplot2)
    
    cbPalette <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    textsize=50
    
    
    density_plot<-ggplot() +
      #<-ggplot(DenDat, aes(x=DenDat$`Normalized count`)) +
      geom_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),bw=bw,size=2,position="identity", alpha=0.6)+
      scale_colour_manual( values = cbPalette)+
      labs(x = "Normalized count",y='density',color='Normalization methods')+
      ggtitle(filenn)+
      theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1.5,"line"))
    
    
    if(typp=='2i'){
      gg_list_2i[[qq2]]<-density_plot
      
    }else if(typp=='S'){
      gg_list_S[[qq2]]<-density_plot
    }

  }
}
names(DenDat_list_2i)<-Fig3Gene
names(DenDat_list_S)<-Fig3Gene


#####begin density plots######

library(gridExtra)
library(ggpubr)
source("E:/RNAseqProject/Density_fun.r")

textsize<-10
bw<-10
gg_list_2i<-lapply(DenDat_list_2i,qwerfun,CAPTION='',linesize=0.5)
gg_list_S<-lapply(DenDat_list_S,qwerfun,CAPTION='',linesize=0.5)
names(gg_list_2i)<-Fig3Gene
names(gg_list_S)<-Fig3Gene
source("E:/RNAseqProject/MANY_SAVE_PATH.R")



gg_list_2i$Pou5f1


###2i whole: Fig S12####
qq<-do.call('ggarrange',c(gg_list_2i,ncol=3,nrow=3,common.legend = TRUE, legend="bottom"))
#dev.off()
ggsave(FIGURE_SUP_PATH_fun('/SF_2i_default_divmeanbeta.pdf'),width = 8, height =8,units='in')

###serum whole: Fig S13####
qq<-do.call('ggarrange',c(gg_list_S,ncol=3,nrow=3,common.legend = TRUE, legend="bottom"))
#dev.off()
ggsave(FIGURE_SUP_PATH_fun('/SF_S_default_divmeanbeta.pdf'),width = 8, height =8,units='in')


library(gridExtra)

qwerfun_ind(DenDat_list_2i$Gli2)

###Stag3#####
textsize<-12
bw<-10
source(file="E:/RNAseqProject/MANY_SAVE_PATH.R")
source("E:/RNAseqProject/Density_fun.r")
qq<-qwerfun_ind(DenDat_list_2i$Stag3,CAPTION='Data from Grün et al study (2i medium)')
qq

pdf("E:/RNAseqProject/Illustrator_bayNorm/newFIG/Grun_Stag3_default_divmeanbeta.pdf",width=6,height=5.3)
qq
dev.off()
