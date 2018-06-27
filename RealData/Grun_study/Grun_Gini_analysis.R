load("Grun_2i_norms.RData")
load("Grun_serum_norms.RData")
load("smFISH_norm_load.RData")

SCnorm_2i<-SCnorm_2i_ori@metadata$NormalizedData
SCnorm_serum<-SCnorm_serum_ori@metadata$NormalizedData


Method_2i<-c(rep("smFISH",9),rep("bayNorm",9),rep("SCnorm",9),rep("Scaling",9),rep("SAVER",9),rep('scImpute',9),rep('MAGIC',9),rep('DCA',9))
Method_serum<-c(rep("smFISH",9),rep("bayNorm",9),rep("SCnorm",9),rep("Scaling",9),rep("SAVER",9),rep('scImpute',9),rep('MAGIC',9),rep('DCA',9))

load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_smFISH_meanBETA/Bootstrapping_Gini_V3tr_default_divmeanbeta.RData")




#scale normalized data by mean BETA, which was used in bayNorm. For a fair comparison.
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



library(ineq)
gini_fish_2i<-unlist(lapply(Fish_2iplot_list,Gini))[Fig3Gene]
gini_bay_2i<-apply(bayNorm_SC_2i$Bay_array[Fig3Gene,,],1,Gini)
gini_scnorm_2i<-apply(SCnorm_2i[Fig3Gene,],1,Gini)
CountsUMI_SC_2i_beta<-t(t(CountsUMI_SC_2i)/bayNorm_SC_2i$BETA)
gini_rawbeta_2i<-apply(CountsUMI_SC_2i_beta[Fig3Gene,],1,Gini)
gini_saver_2i<-apply(SAVER_SC_2i_array[Fig3Gene,,],1,Gini)
gini_scimpute_2i<-apply(scImpute_SC_2i[Fig3Gene,],1,Gini)
gini_magic_2i<-apply(MAGIC_2i[Fig3Gene,],1,Gini)
gini_dca_2i<-apply(DCA_2i[Fig3Gene,],1,Gini)


gini_fish_serum<-unlist(lapply(Fish_serumplot_list,Gini))[Fig3Gene]
gini_bay_serum<-apply(bayNorm_SC_serum$Bay_array[Fig3Gene,,],1,Gini)
gini_scnorm_serum<-apply(SCnorm_serum[Fig3Gene,],1,Gini)
CountsUMI_SC_serum_beta<-t(t(CountsUMI_SC_serum)/bayNorm_SC_serum$BETA)
gini_rawbeta_serum<-apply(CountsUMI_SC_serum_beta[Fig3Gene,],1,Gini)
gini_saver_serum<-apply(SAVER_SC_serum_array[Fig3Gene,,],1,Gini)
gini_scimpute_serum<-apply(scImpute_SC_serum[Fig3Gene,],1,Gini)
gini_magic_serum<-apply(MAGIC_serum[Fig3Gene,],1,Gini)
gini_dca_serum<-apply(DCA_serum[Fig3Gene,],1,Gini)



#####smFISH#####
ERROR_2i_smFISH<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_smFISH,2,var))[Fig3Gene]
ERROR_serum_smFISH<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_smFISH,2,var))[Fig3Gene]
MEAN_2i_smFISH<-colMeans(STORE_BOOTGINI_list$boot_2imat_smFISH)[Fig3Gene]
MEAN_serum_smFISH<-colMeans(STORE_BOOTGINI_list$boot_serummat_smFISH)[Fig3Gene]

#####bayNorm#####
ERROR_2i_bayNorm<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_bayNorm,2,var))[Fig3Gene]
ERROR_serum_bayNorm<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_bayNorm,2,var))[Fig3Gene]
MEAN_2i_bayNorm<-colMeans(STORE_BOOTGINI_list$boot_2imat_bayNorm)[Fig3Gene]
MEAN_serum_bayNorm<-colMeans(STORE_BOOTGINI_list$boot_serummat_bayNorm)[Fig3Gene]

ERROR_2i_scnorm<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_scnorm,2,var))[Fig3Gene]
ERROR_serum_scnorm<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_scnorm,2,var))[Fig3Gene]
MEAN_2i_scnorm<-colMeans(STORE_BOOTGINI_list$boot_2imat_scnorm)[Fig3Gene]
MEAN_serum_scnorm<-colMeans(STORE_BOOTGINI_list$boot_serummat_scnorm)[Fig3Gene]

ERROR_2i_scran<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_scran,2,var))[Fig3Gene]
ERROR_serum_scran<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_scran,2,var))[Fig3Gene]
MEAN_2i_scran<-colMeans(STORE_BOOTGINI_list$boot_2imat_scran)[Fig3Gene]
MEAN_serum_scran<-colMeans(STORE_BOOTGINI_list$boot_serummat_scran)[Fig3Gene]

ERROR_2i_RPM<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_RPM,2,var))[Fig3Gene]
ERROR_serum_RPM<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_RPM,2,var))[Fig3Gene]
MEAN_2i_RPM<-colMeans(STORE_BOOTGINI_list$boot_2imat_RPM)[Fig3Gene]
MEAN_serum_RPM<-colMeans(STORE_BOOTGINI_list$boot_serummat_RPM)[Fig3Gene]

ERROR_2i_TMM<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_TMM,2,var))[Fig3Gene]
ERROR_serum_TMM<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_TMM,2,var))[Fig3Gene]
MEAN_2i_TMM<-colMeans(STORE_BOOTGINI_list$boot_2imat_TMM)[Fig3Gene]
MEAN_serum_TMM<-colMeans(STORE_BOOTGINI_list$boot_serummat_TMM)[Fig3Gene]


ERROR_2i_DESeq<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_DESeq,2,var))[Fig3Gene]
ERROR_serum_DESeq<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_DESeq,2,var))[Fig3Gene]
MEAN_2i_DESeq<-colMeans(STORE_BOOTGINI_list$boot_2imat_DESeq)[Fig3Gene]
MEAN_serum_DESeq<-colMeans(STORE_BOOTGINI_list$boot_serummat_DESeq)[Fig3Gene]

ERROR_2i_SAVER<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_SAVER,2,var))[Fig3Gene]
ERROR_serum_SAVER<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_SAVER,2,var))[Fig3Gene]
MEAN_2i_SAVER<-colMeans(STORE_BOOTGINI_list$boot_2imat_SAVER)[Fig3Gene]
MEAN_serum_SAVER<-colMeans(STORE_BOOTGINI_list$boot_serummat_SAVER)[Fig3Gene]



ERROR_2i_RAW_BETA<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_RAW_BETA,2,var))[Fig3Gene]
ERROR_serum_RAW_BETA<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_RAW_BETA,2,var))[Fig3Gene]
MEAN_2i_RAW_BETA<-colMeans(STORE_BOOTGINI_list$boot_2imat_RAW_BETA)[Fig3Gene]
MEAN_serum_RAW_BETA<-colMeans(STORE_BOOTGINI_list$boot_serummat_RAW_BETA)[Fig3Gene]


ERROR_2i_scImpute<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_scImpute,2,var))[Fig3Gene]
ERROR_serum_scImpute<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_scImpute,2,var))[Fig3Gene]
MEAN_2i_scImpute<-colMeans(STORE_BOOTGINI_list$boot_2imat_scImpute)[Fig3Gene]
MEAN_serum_scImpute<-colMeans(STORE_BOOTGINI_list$boot_serummat_scImpute)[Fig3Gene]


ERROR_2i_MAGIC<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_MAGIC,2,var))[Fig3Gene]
ERROR_serum_MAGIC<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_MAGIC,2,var))[Fig3Gene]
MEAN_2i_MAGIC<-colMeans(STORE_BOOTGINI_list$boot_2imat_MAGIC)[Fig3Gene]
MEAN_serum_MAGIC<-colMeans(STORE_BOOTGINI_list$boot_serummat_MAGIC)[Fig3Gene]

ERROR_2i_DCA<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_2imat_DCA,2,var))[Fig3Gene]
ERROR_serum_DCA<-1.96*sqrt(apply(STORE_BOOTGINI_list$boot_serummat_DCA,2,var))[Fig3Gene]
MEAN_2i_DCA<-colMeans(STORE_BOOTGINI_list$boot_2imat_DCA)[Fig3Gene]
MEAN_serum_DCA<-colMeans(STORE_BOOTGINI_list$boot_serummat_DCA)[Fig3Gene]

Gini_2i<-c(gini_fish_2i[Fig3Gene],gini_bay_2i[Fig3Gene],gini_scnorm_2i[Fig3Gene],gini_rawbeta_2i[Fig3Gene],gini_saver_2i[Fig3Gene],gini_scimpute_2i[Fig3Gene],gini_magic_2i[Fig3Gene],gini_dca_2i[Fig3Gene])
ERROR_2i<-c(ERROR_2i_smFISH[Fig3Gene],ERROR_2i_bayNorm[Fig3Gene],ERROR_2i_scnorm[Fig3Gene],ERROR_2i_RAW_BETA[Fig3Gene],ERROR_2i_SAVER[Fig3Gene],ERROR_2i_scImpute[Fig3Gene],ERROR_2i_MAGIC[Fig3Gene],ERROR_2i_DCA[Fig3Gene])

MEAN_2i<-c(MEAN_2i_smFISH[Fig3Gene],MEAN_2i_bayNorm[Fig3Gene],MEAN_2i_scnorm[Fig3Gene],MEAN_2i_RAW_BETA[Fig3Gene],MEAN_2i_SAVER[Fig3Gene],MEAN_2i_scImpute[Fig3Gene],MEAN_2i_MAGIC[Fig3Gene],MEAN_2i_DCA[Fig3Gene])


BAR_DAT_2i<-data.frame(Gini=Gini_2i,Method=Method_2i,Genename=rep(Fig3Gene,length(unique(Method_2i))),ERROR=ERROR_2i,MEAN=MEAN_2i)



Gini_serum<-c(gini_fish_serum[Fig3Gene],gini_bay_serum[Fig3Gene],gini_scnorm_serum[Fig3Gene],gini_rawbeta_serum[Fig3Gene],gini_saver_serum[Fig3Gene],gini_scimpute_serum[Fig3Gene],gini_magic_serum[Fig3Gene],gini_dca_serum[Fig3Gene])

ERROR_serum<-c(ERROR_serum_smFISH[Fig3Gene],ERROR_serum_bayNorm[Fig3Gene],ERROR_serum_scnorm[Fig3Gene],ERROR_serum_RAW_BETA[Fig3Gene],ERROR_serum_SAVER[Fig3Gene],ERROR_serum_scImpute[Fig3Gene],ERROR_serum_MAGIC[Fig3Gene],ERROR_serum_DCA[Fig3Gene])

MEAN_serum<-c(MEAN_serum_smFISH[Fig3Gene],MEAN_serum_bayNorm[Fig3Gene],MEAN_serum_scnorm[Fig3Gene],MEAN_serum_RAW_BETA[Fig3Gene],MEAN_serum_SAVER[Fig3Gene],MEAN_serum_scImpute[Fig3Gene],MEAN_serum_MAGIC[Fig3Gene],MEAN_serum_DCA[Fig3Gene])


BAR_DAT_S<-data.frame(Gini=Gini_serum,Method=Method_serum,Genename=rep(Fig3Gene,length(unique(Method_serum))),ERROR=ERROR_serum,MEAN=MEAN_serum)



#BEGIN PLOT!!!!!!!!!!!!!!!!!
#choose a medium
BAR_DAT_input<-BAR_DAT_2i
CAPTION<-'Single cell data based on 2i medium'

BAR_DAT_input<-BAR_DAT_S
CAPTION<-'Single cell data based on serum medium'

BAR_DAT_input[,1]<-as.numeric(as.character(BAR_DAT_input[,1]))
BAR_DAT_input[,3]<-factor(BAR_DAT_input[,3],levels=unique(BAR_DAT_input[,3]))
BAR_DAT_input[,2]<-factor(BAR_DAT_input[,2],levels=unique(BAR_DAT_input[,2]))



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


library(ggplot2)
textsize<-12
BAR_out<-ggplot(data=BAR_DAT_input, aes(x=BAR_DAT_input[,3], y=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN, fill=BAR_DAT_input[,2])) +
  geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
  geom_errorbar(aes(ymin=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN-BAR_DAT_input$ERROR, ymax=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN+BAR_DAT_input$ERROR),size=.5,width=.4,position=position_dodge(.9))+
  #geom_text(aes(label=round(2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN,2)), vjust=1.6, color="black", position = position_dodge(0.9), size=2.5)+
  labs(x = "Gene",y='Gini',fill='Methods',caption=CAPTION)+
  #ggtitle(CAPTION) +
  #scale_fill_brewer(palette='Paired)+
  scale_fill_manual(values=cbbPalette)+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position = 'top')
#dev.off()
#legend.position =c(0.9,0.85)

BAR_out

#pepare for Figure 2####
BAR_DAT_2i_GINI<-BAR_DAT_2i
BAR_DAT_S_GINI<-BAR_DAT_S
save(BAR_DAT_2i_GINI,BAR_DAT_S_GINI,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_smFISH_meanBETA/GINI_data_default_divmeanbeta.RData")




