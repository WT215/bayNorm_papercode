source("~/SIM_FigS15_functions.R")


source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")
#SIM1: B_01_01_D_02_0#############
load("E:/RNAseqProject/SIMULATION/SIM_1/SIM_1.RData")



INPUT_LIST_g1<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 1`,SCnorm=scnorm_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),Scaling=RB_out[,seq(1,100)],SAVER=saver_array1/mean(bayNorm_out$BETA[[1]]),scImpute=scImpute_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),MAGIC=MAGIC_SIM_1[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),DCA=DCA_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]))

INPUT_LIST_g2<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 2`,SCnorm=scnorm_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),Scaling=RB_out[,-seq(1,100)],SAVER=saver_array2/mean(bayNorm_out$BETA[[2]]),scImpute=scImpute_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),MAGIC=MAGIC_SIM_1[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),DCA=DCA_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]))


RE_1_g1<-SIM_ana_fun(INPUT_LIST_g1)
RE_1_g2<-SIM_ana_fun(INPUT_LIST_g2)
temp_list<-list(RE_1_g1,RE_1_g2)


TRUE_MU_S1<-TrueMU
DROP_S1<-DROP
TRUE_DAT_S1<-SCE@assays$data$TrueCounts[-DROP_S1,]

TRUE_MU_S1_1<-rowMeans(TRUE_DAT_S1[,CONDITION==1])
TRUE_MU_S1_2<-rowMeans(TRUE_DAT_S1[,CONDITION==2])
TRUE_CV_S1_1<-CV_fun(TRUE_DAT_S1[,CONDITION==1])
TRUE_CV_S1_2<-CV_fun(TRUE_DAT_S1[,CONDITION==2])
TRUE_Gini_S1_1<-apply(TRUE_DAT_S1[,CONDITION==1],1,Gini)
TRUE_Gini_S1_2<-apply(TRUE_DAT_S1[,CONDITION==2],1,Gini)



textsize<-6
lwd=0.05
outlier.size<-0.001


mean_re_S1<-mean_data_fun(temp_list=temp_list,TRUE_MU =TRUE_MU_S1,textsize = textsize,lwd=lwd,outlier.size=outlier.size)

mean_re_S1$MEAN_PLOT


cv_re_S1<-cv_data_fun(temp_list=temp_list,textsize = textsize,TRUE_CV_g1=TRUE_CV_S1_1,TRUE_CV_g2=TRUE_CV_S1_2,lwd=lwd,outlier.size=outlier.size)
gini_re_S1<-gini_data_fun(temp_list=temp_list,textsize = textsize,TRUE_Gini_g1=TRUE_Gini_S1_1,TRUE_Gini_g2=TRUE_Gini_S1_2,lwd=lwd,outlier.size=outlier.size)

#SIM2: B_005_01_D_02_0#############
load("E:/RNAseqProject/SIMULATION/SIM_005_01/SIM_005_01.RData")



INPUT_LIST_g1<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 1`,SCnorm=scnorm_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),Scaling=RB_out[,seq(1,100)],SAVER=saver_array1/mean(bayNorm_out$BETA[[1]]),scImpute=scImpute_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),MAGIC=MAGIC_SIM_005_01[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),DCA=DCA_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]))

INPUT_LIST_g2<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 2`,SCnorm=scnorm_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),Scaling=RB_out[,-seq(1,100)],SAVER=saver_array2/mean(bayNorm_out$BETA[[2]]),scImpute=scImpute_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),MAGIC=MAGIC_SIM_005_01[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),DCA=DCA_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]))

source("E:/RNAseqProject/SIMULATION/SIM_mean_cv_gini_fun.R")
RE_005_01_g1<-SIM_ana_fun(INPUT_LIST_g1)
RE_005_01_g2<-SIM_ana_fun(INPUT_LIST_g2)
temp_list<-list(RE_005_01_g1,RE_005_01_g2)


TRUE_MU_S2<-TrueMU
DROP_S2<-DROP

TRUE_DAT_S2<-SCE@assays$data$TrueCounts[-DROP_S2,]

TRUE_MU_S2_1<-rowMeans(TRUE_DAT_S2[,CONDITION==1])
TRUE_MU_S2_2<-rowMeans(TRUE_DAT_S2[,CONDITION==2])
TRUE_CV_S2_1<-CV_fun(TRUE_DAT_S2[,CONDITION==1])
TRUE_CV_S2_2<-CV_fun(TRUE_DAT_S2[,CONDITION==2])
TRUE_Gini_S2_1<-apply(TRUE_DAT_S2[,CONDITION==1],1,Gini)
TRUE_Gini_S2_2<-apply(TRUE_DAT_S2[,CONDITION==2],1,Gini)



mean_re_S2<-mean_data_fun(temp_list=temp_list,TRUE_MU =TRUE_MU_S2,textsize = textsize,lwd=lwd,outlier.size=outlier.size)


cv_re_S2<-cv_data_fun(temp_list=temp_list,textsize = textsize,TRUE_CV_g1=TRUE_CV_S2_1,TRUE_CV_g2=TRUE_CV_S2_2,lwd=lwd,outlier.size=outlier.size)
gini_re_S2<-gini_data_fun(temp_list=temp_list,textsize = textsize,TRUE_Gini_g1=TRUE_Gini_S2_1,TRUE_Gini_g2=TRUE_Gini_S2_2,lwd=lwd,outlier.size=outlier.size)





#SIM3: B_005_005_D_02_0#############
load("E:/RNAseqProject/SIMULATION/SIM_005_005/SIM_005_005.RData")


INPUT_LIST_g1<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 1`,SCnorm=scnorm_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),Scaling=RB_out[,seq(1,100)],SAVER=saver_array1/mean(bayNorm_out$BETA[[1]]),scImpute=scImpute_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),MAGIC=MAGIC_SIM_005_005[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),DCA=DCA_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]))

INPUT_LIST_g2<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 2`,SCnorm=scnorm_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),Scaling=RB_out[,-seq(1,100)],SAVER=saver_array2/mean(bayNorm_out$BETA[[2]]),scImpute=scImpute_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),MAGIC=MAGIC_SIM_005_005[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),DCA=DCA_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]))

source("E:/RNAseqProject/SIMULATION/SIM_mean_cv_gini_fun.R")
RE_005_005_g1<-SIM_ana_fun(INPUT_LIST_g1)
RE_005_005_g2<-SIM_ana_fun(INPUT_LIST_g2)
temp_list<-list(RE_005_005_g1,RE_005_005_g2)


TRUE_MU_S3<-TrueMU
DROP_S3<-DROP

TRUE_DAT_S3<-SCE@assays$data$TrueCounts[-DROP_S3,]

TRUE_MU_S3_1<-rowMeans(TRUE_DAT_S3[,CONDITION==1])
TRUE_MU_S3_2<-rowMeans(TRUE_DAT_S3[,CONDITION==2])
TRUE_CV_S3_1<-CV_fun(TRUE_DAT_S3[,CONDITION==1])
TRUE_CV_S3_2<-CV_fun(TRUE_DAT_S3[,CONDITION==2])
TRUE_Gini_S3_1<-apply(TRUE_DAT_S3[,CONDITION==1],1,Gini)
TRUE_Gini_S3_2<-apply(TRUE_DAT_S3[,CONDITION==2],1,Gini)


mean_re_S3<-mean_data_fun(temp_list=temp_list,TRUE_MU =TRUE_MU_S3,textsize = textsize,lwd=lwd,outlier.size=outlier.size)


cv_re_S3<-cv_data_fun(temp_list=temp_list,textsize = textsize,TRUE_CV_g1=TRUE_CV_S3_1,TRUE_CV_g2=TRUE_CV_S3_2,lwd=lwd,outlier.size=outlier.size)
gini_re_S3<-gini_data_fun(temp_list=temp_list,textsize = textsize,TRUE_Gini_g1=TRUE_Gini_S3_1,TRUE_Gini_g2=TRUE_Gini_S3_2,lwd=lwd,outlier.size=outlier.size)


mean_re_S3$MEAN_PLOT
cv_re_S3$CV_PLOT
gini_re_S3$Gini_PLOT




#SIM4: B_01_005_D_02_0#############
load("E:/RNAseqProject/SIMULATION/SIM_01_005/SIM_01_005.RData")


INPUT_LIST_g1<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 1`,SCnorm=scnorm_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),Scaling=RB_out[,seq(1,100)],SAVER=saver_array1/mean(bayNorm_out$BETA[[1]]),scImpute=scImpute_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),MAGIC=MAGIC_SIM_01_005[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]),DCA=DCA_out[,seq(1,100)]/mean(bayNorm_out$BETA[[1]]))

INPUT_LIST_g2<-list(bayNorm=bayNorm_out$Bay_array_list$`Group 2`,SCnorm=scnorm_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),Scaling=RB_out[,-seq(1,100)],SAVER=saver_array2/mean(bayNorm_out$BETA[[2]]),scImpute=scImpute_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),MAGIC=MAGIC_SIM_01_005[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]),DCA=DCA_out[,-seq(1,100)]/mean(bayNorm_out$BETA[[2]]))


RE_01_005_g1<-SIM_ana_fun(INPUT_LIST_g1)
RE_01_005_g2<-SIM_ana_fun(INPUT_LIST_g2)
temp_list<-list(RE_01_005_g1,RE_01_005_g2)


TRUE_MU_S4<-TrueMU
DROP_S4<-DROP

TRUE_DAT_S4<-SCE@assays$data$TrueCounts[-DROP_S4,]

TRUE_MU_S4_1<-rowMeans(TRUE_DAT_S4[,CONDITION==1])
TRUE_MU_S4_2<-rowMeans(TRUE_DAT_S4[,CONDITION==2])
TRUE_CV_S4_1<-CV_fun(TRUE_DAT_S4[,CONDITION==1])
TRUE_CV_S4_2<-CV_fun(TRUE_DAT_S4[,CONDITION==2])
TRUE_Gini_S4_1<-apply(TRUE_DAT_S4[,CONDITION==1],1,Gini)
TRUE_Gini_S4_2<-apply(TRUE_DAT_S4[,CONDITION==2],1,Gini)


mean_re_S4<-mean_data_fun(temp_list=temp_list,TRUE_MU =TRUE_MU_S4,textsize = textsize,lwd=lwd,outlier.size=outlier.size)


cv_re_S4<-cv_data_fun(temp_list=temp_list,textsize = textsize,TRUE_CV_g1=TRUE_CV_S4_1,TRUE_CV_g2=TRUE_CV_S4_2,lwd=lwd,outlier.size=outlier.size)
gini_re_S4<-gini_data_fun(temp_list=temp_list,textsize = textsize,TRUE_Gini_g1=TRUE_Gini_S4_1,TRUE_Gini_g2=TRUE_Gini_S4_2,lwd=lwd,outlier.size=outlier.size)




###Begin plot for Fig S15####

library(gridExtra)
library(ggpubr)
library(cowplot)

qq<-ggarrange(
    mean_re_S1$MEAN_PLOT + theme(legend.position='none'),
    cv_re_S1$CV_PLOT + theme(legend.position="none"),
    gini_re_S1$Gini_PLOT+ theme(legend.position="none"),
    
    mean_re_S2$MEAN_PLOT + theme(legend.position="none"),
    cv_re_S2$CV_PLOT + theme(legend.position="none"),
    gini_re_S2$Gini_PLOT+ theme(legend.position="none"), 
    
    mean_re_S4$MEAN_PLOT + theme(legend.position="none"),
    cv_re_S4$CV_PLOT + theme(legend.position="none"),
    gini_re_S4$Gini_PLOT+ theme(legend.position="none"), 
    
    mean_re_S3$MEAN_PLOT + theme(legend.position="none"),
    cv_re_S3$CV_PLOT + theme(legend.position="none"),
    gini_re_S3$Gini_PLOT+ theme(legend.position="none"), 

    ncol=3,nrow=4,common.legend = TRUE, legend="top")
qq

ggsave(filename ="E:/RNAseqProject/Illustrator_bayNorm/newFIG/SUP_SIM_meancvgini.pdf",plot = qq,width=11,height=8)
