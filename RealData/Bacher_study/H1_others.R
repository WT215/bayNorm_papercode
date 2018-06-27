load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/meanBay/H1_BAY_esf_01.RData")
load("E:/RNAseqProject/Bacher__SCnorm_2016/RAW_INITIATE.RData")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")



H1_data_comb<-cbind(H1_p24,H1_p96)[whichg_H1,]
dim(H1_data_comb)



library(bayNorm)


aBAY_H1<-bayNorm_sup(Data=H1_data_comb,Conditions=CONDITION_H1,UMI_sffl = c(20,10),S=10,PRIORS = mBAY_H1$PRIORS_LIST,BETA_vec = do.call(c,mBAY_H1$BETA))


library(foreach)

M_array_mat<-foreach(j=1:10,.combine=cbind)%do%{
    print(j)
    qq<-SCnorm_runMAST(Data=cbind(aBAY_H1$Bay_array_list$`Group 1`[,,j],aBAY_H1$Bay_array_list$`Group 2`[,,j]), NumCells=c(92,92))$adjpval
    return(qq)
}

library(ROTS)
R_array_mat<-foreach(j=1:10,.combine=cbind)%do%{
    print(j)
    qq<- ROTS(data=cbind(aBAY_H1$Bay_array_list$`Group 1`[,,j],aBAY_H1$Bay_array_list$`Group 2`[,,j]),groups=CONDITION_H1,B=100,log = F)$FDR
    names(qq)<-rownames(aBAY_H1$Bay_array_list$`Group 1`)
    return(qq)
}
#for each sample run wilcox, then take median
W_array<-wilcox_fun(cells=aBAY_H1$Bay_array_list$`Group 1`,ctrls=aBAY_H1$Bay_array_list$`Group 2`)

save(aBAY_H1,M_array_mat,R_array_mat,W_array,file="E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/aBAY_H1.RData")


######MODE############

#modeBAY_H1<-bayNorm_sup(Data=H1_data_comb,Conditions=CONDITION_H1,UMI_sffl = c(20,10),S=10,mode_version = T,PRIORS = mBAY_H1$PRIORS_LIST,BETA_vec = do.call(c,mBAY_H1$BETA))

M_mode<-SCnorm_runMAST(Data=do.call(cbind,modeBAY_H1$Bay_mat_list), NumCells=c(92,92))$adjpval

R_mode<-ROTS(data=do.call(cbind,modeBAY_H1$Bay_mat_list),groups=CONDITION_H1,B=100,log = F)$FDR
names(R_mode)<-rownames(modeBAY_H1$Bay_mat_list$`Group 1`)

W_mode<-wilcox_fun(cells=modeBAY_H1$Bay_mat_list$`Group 1`,ctrls=modeBAY_H1$Bay_mat_list$`Group 2`)


save(modeBAY_H1,M_mode,R_mode,W_mode,file="E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/modeBAY_H1.RData")

#####MEAN###########
M_mean<-M_mBAY_H1$adjpval

R_mean<-ROTS(data=do.call(cbind,mBAY_H1$Bay_mat_list),groups=CONDITION_H1,B=100,log = F)$FDR
names(R_mean)<-rownames(aBAY_H1$Bay_array_list$`Group 1`)

W_mean<-wilcox_fun(cells=mBAY_H1$Bay_mat_list$`Group 1`,ctrls=mBAY_H1$Bay_mat_list$`Group 2`)

save(mBAY_H1,M_mean,R_mean,W_mean,file="E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/meanBAY_H1.RData")


#######begin plot#########
load("E:/RNAseqProject/Bacher__SCnorm_2016/RAW_INITIATE.RData")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/aBAY_H1.RData")
load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/modeBAY_H1.RData")
load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/H1_others/meanBAY_H1.RData")

Inputdat<-H1_data_comb
MedExp <- log(apply(Inputdat, 1, function(x) median(x[x != 0])))
grpnum <- 6
splitby <- sort(MedExp)
grps <- length(splitby)/grpnum
sreg <- split(splitby, ceiling(seq_along(splitby)/grps))

Input_re_list<-list(M_mean,R_mean,W_mean,M_mode,R_mode,W_mode,apply(M_array_mat,1,median),apply(R_array_mat,1,median),W_array)

library(foreach)
method_vec<-c('M_mean','R_mean','W_mean','M_mode','R_mode','W_mode','M_array','R_array','W_array')
type_vec<-c(rep('mean',3),rep('mode',3),rep('array',3))
DEmethod_vec<-rep(c('MAST','ROTS','Wilcoxon'),3)


INPUT_LIST<-Input_re_list
names(INPUT_LIST)<-method_vec

Gene_exp_gr<-seq(1,grpnum)

library(foreach)
thres<-0.05

BAR_DAT<-foreach(i=1:length(sreg),.combine=rbind)%:%
    foreach(j=1:length(INPUT_LIST),.combine=rbind)%do%{
        qq<-length(intersect(names(which(INPUT_LIST[[j]]<thres)),names(sreg[[i]])))
        qq2<-c(qq,names(INPUT_LIST)[j],i)
        return(qq2)
        
    }
    


BAR_DAT<-as.data.frame(BAR_DAT)
colnames(BAR_DAT)<-c('Number of detected DE genes','Normalization methods','Gene expression group')
BAR_DAT[,1]<-as.numeric(as.character(BAR_DAT[,1]))
BAR_DAT[,2]<-factor(BAR_DAT[,2],levels=unique(BAR_DAT[,2]))
BAR_DAT[,3]<-factor(BAR_DAT[,3],levels=unique(BAR_DAT[,3]))


library(ggplot2)
textsize<-12
DE_H1<-ggplot(data=BAR_DAT, aes(x=BAR_DAT[,2], y=BAR_DAT[,1], fill=BAR_DAT[,3])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=BAR_DAT[,1]), vjust=0, color="black", position = position_dodge(0.9), size=2)+
    labs(fill = "Gene expression group",y='Number of detected DE genes',x='DE method_type of output from bayNorm')+ggtitle("") +
    scale_fill_brewer(palette="Paired")+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
#dev.off()
DE_H1

source(file="E:/RNAseqProject/MANY_SAVE_PATH.R")
ggsave(filename=FIGURE_SUP_PATH_fun('/SUP_H1_others.pdf'),plot=DE_H1,width=8.2,height=5.5)



library(cowplot)
library(gridExtra)
library(ggpubr)
qqq<-ggarrange(DE_H1,BAR_OUT, BAR_OUT_simH1,ncol=1, nrow=3, common.legend = TRUE, legend="top")
qqq

ggsave(filename=FIGURE_SUP_PATH_fun('/SUP_bayDE.pdf'),plot=qqq,width=8.2,height=9.5)
# qqq<-plot_grid(DE_H1+ theme(legend.position="none"),
#                BAR_OUT+ theme(legend.position="none"),
#                nrow=2,ncol=1)

