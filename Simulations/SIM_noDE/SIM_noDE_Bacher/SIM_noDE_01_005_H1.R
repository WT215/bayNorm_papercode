source("E:/RNAseqProject/QQsim_v2_SingleCellExperiment.R")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

load("E:/RNAseqProject/Bacher__SCnorm_2016/DIY_sim/H1p24_bay_sim_allgene.RData")


Sim_List_Input<-H1p24_bay_sim

Sim_List_Input$SCE<-NULL
Sim_List_Input$Est_params@counts.norm.TC<-matrix()
library(bayNorm)


#######Begin simulation#######


Est_params<-H1p24_bay_sim$Est_params
Est_params@nGroups<-2
Est_params@groupCells=c(100,100)
Est_params@nGenes<-10000
Est_params@MeanBeta<-c(0.1,0.05)
Est_params@de.prob<-c(0,0)
Est_params@de.downProb<-c(0,0)
Est_params@de.facLoc<-c(1,1)
Est_params@de.facScale<-c(0.5,0.5)

Est_params@out.prob

Est_params@beta.scale

Est_params@out.prob
Est_params@mean.shape

Est_params@Beta<-numeric()
length(Est_params@Beta)


library(splatter)


SCE<-QQinitiate(Est_params)
SCE<-QQSimGeneMeans(SCE,Est_params)
#Group
SCE<-QQSimGroupDE(SCE,Est_params)
SCE<-QQSimGroupCellMeans(SCE, Est_params,Effect=F)

SCE<-QQSimBETA(SCE, Est_params)
SCE<-QQSimBCVMeans(SCE,Est_params)
SCE<-QQSimBinomial(SCE,Est_params)


dim(counts(SCE))
CONDITION<-SCE@colData@listData$GroupInd
table(CONDITION)


DROP<-which(rowSums(counts(SCE))==0)

Est_params@mean.shape
Est_params@mean.rate
Est_params@out.prob
Est_params@out.facLoc
Est_params@out.facScale
Est_params@bcv.common
Est_params@bcv.df
Est_params@beta.loc
Est_params@beta.scale

#########True Parameters#####
TrueBeta1<-colData(SCE)$Beta[CONDITION==1]*Est_params@MeanBeta[1]
TrueBeta2<-colData(SCE)$Beta[CONDITION==2]*Est_params@MeanBeta[2]
summary(TrueBeta1)
summary(TrueBeta2)




TrueSIZES1<-1/(rowMeans(SCE@assays@.xData$data$BCV[,CONDITION==1])[-DROP])^2
TrueSIZES2<-1/(rowMeans(SCE@assays@.xData$data$BCV[,CONDITION==2])[-DROP])^2

TrueMU<-rowData(SCE)$BaseGeneMean[-DROP]


#begin bayNorm######
RAW_DAT<-counts(SCE)[-DROP,]
library(scran)
BETA_SCRAN_1<-computeSumFactors(RAW_DAT[,CONDITION==1])
BETA_SCRAN_1<-BETA_SCRAN_1/mean(BETA_SCRAN_1)*Est_params@MeanBeta[1]
BETA_SCRAN_2<-computeSumFactors(RAW_DAT[,CONDITION==2])
BETA_SCRAN_2<-BETA_SCRAN_2/mean(BETA_SCRAN_2)*Est_params@MeanBeta[2]

plot(BETA_SCRAN_1,TrueBeta1,pch=16)
abline(0,1)


library(bayNorm)
#bayNorm_out<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,Prior_type = 'GG',S=10)

mbayNorm_out<-bayNorm(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,mean_version = T,S=1000,Prior_type = 'GG')


gr<-table(CONDITION)
M_mbayNorm<-SCnorm_runMAST(Data=cbind(mbayNorm_out$Bay_mat_list$`Group 1`,mbayNorm_out$Bay_mat_list$`Group 2`),NumCells = gr)
T_mbayNorm<-wilcox_fun(cells=mbayNorm_out$Bay_mat_list$`Group 1`,ctrls=mbayNorm_out$Bay_mat_list$`Group 2`)

R_mbayNorm<-ROTS(data=cbind(mbayNorm_out$Bay_mat_list$`Group 1`,mbayNorm_out$Bay_mat_list$`Group 2`),groups=CONDITION,B=100,log=F)
names(R_mbayNorm)<-rownames(RAW_DAT)






#mode 
modebayNorm_out<-bayNorm_sup(Data=RAW_DAT,BETA_vec=do.call(c,mbayNorm_out$BETA),Conditions=CONDITION,mode_version = T,S=1000,PRIORS=mbayNorm_out$PRIORS_LIST)

modebayNorm_out<-bayNorm_sup(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,mode_version = T,S=1000,PRIORS =mbayNorm_out$PRIORS_LIST )
M_modebayNorm<-SCnorm_runMAST(Data=cbind(modebayNorm_out$Bay_mat_list$`Group 1`,modebayNorm_out$Bay_mat_list$`Group 2`),NumCells = gr)
T_modebayNorm<-wilcox_fun(cells=modebayNorm_out$Bay_mat_list$`Group 1`,ctrls=modebayNorm_out$Bay_mat_list$`Group 2`)
library(ROTS)
R_modebayNorm<-ROTS(data=cbind(modebayNorm_out$Bay_mat_list$`Group 1`,modebayNorm_out$Bay_mat_list$`Group 2`),groups=CONDITION,B=100,log=F)
names(R_modebayNorm)<-rownames(RAW_DAT)

#array
abayNorm_out<-bayNorm_sup(Data=RAW_DAT,BETA_vec=do.call(c,mbayNorm_out$BETA),Conditions=CONDITION,S=10,PRIORS=mbayNorm_out$PRIORS_LIST)

abayNorm_out<-bayNorm_sup(Data<-RAW_DAT,BETA_vec<-c(BETA_SCRAN_1,BETA_SCRAN_2),Conditions=CONDITION,S=10,PRIORS =mbayNorm_out$PRIORS_LIST )
library(abind)


qq<-abind(modebayNorm_out$Bay_array_list,along=2)
M_abayNorm<-SCnorm_runMAST3(Data=qq,NumCells = gr)
T_abayNorm<-wilcox_fun(cells=abayNorm_out$Bay_array_list$`Group 1`,ctrls=abayNorm_out$Bay_array_list$`Group 2`)
library(ROTS)
library(foreach)
R_abayNorm_temp<-foreach(i=1:10,.combine=cbind)%do%{
    qq<-ROTS(data=cbind(modebayNorm_out$Bay_mat_list$`Group 1`,modebayNorm_out$Bay_mat_list$`Group 2`),groups=CONDITION,B=100,log=F)
    names(qq)<-rownames(RAW_DAT)
    return(qq)
}
R_abayNorm<-apply(R_abayNorm_temp,1,median) 

save.image("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")



####SAVER#######
load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")
library(SAVER)
saver_out<-saver(x=RAW_DAT,size.factor=1/c(BETA_SCRAN_1,BETA_SCRAN_2))

gr<-table(CONDITION)
M_saver<-SCnorm_runMAST(Data=saver_out$estimate,NumCells = gr)
#T_saver<-wilcox_fun(cells=saver_out$estimate[,CONDITION==1],ctrls=saver_out$estimate[,CONDITION==2])


save(saver_out,M_saver,file="E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SAVER_noDE_01_005_H1.RData")


#SCnorm
library(SCnorm)
scnorm_out_ori<-SCnorm(Data=RAW_DAT,Conditions=CONDITION,ditherCounts =T)
scnorm_out<-scnorm_out_ori@metadata$NormalizedData
M_scnorm<-SCnorm_runMAST(Data=scnorm_out,NumCells = gr)
T_scnorm<-wilcox_fun(cells=scnorm_out[,CONDITION==1],ctrls=scnorm_out[,CONDITION==2])
save(scnorm_out_ori,scnorm_out,M_scnorm,T_scnorm,file="E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SCnorm_noDE_01_005_H1.RData")


#scImpute
library(scImpute)
scimpute_fun(Data=RAW_DAT,Data_name='DE_noDE_01_005_H1')
scImpute_out<-readRDS("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/DE_noDE_01_005_H1scimpute_count.rds")
gr<-table(CONDITION)
M_scImpute<-SCnorm_runMAST(Data=scImpute_out,NumCells = gr)
T_scImpute<-wilcox_fun(cells=scImpute_out[,which(CONDITION==1)],ctrls=scImpute_out[,which(CONDITION==2)])

save(scImpute_out,M_scImpute,T_scImpute,file="E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/scImpute_noDE_01_005_H1.RData")

#scaling
RB_out<-t(t(RAW_DAT)/c(BETA_SCRAN_1,BETA_SCRAN_2))
gr<-table(CONDITION)
M_RB<-SCnorm_runMAST(Data=RB_out,NumCells = gr)
T_RB<-wilcox_fun(cells=RB_out[,which(CONDITION==1)],ctrls=RB_out[,which(CONDITION==2)])

save(RB_out,M_RB,T_RB,file="E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/RB_noDE_01_005_H1.RData")


##MAGIC#######
load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")

qq<-prcomp(t(MAGIC_out),scale. = T)
plot(qq$x[,1:2],col=CONDITION)

library(Rmagic)
MAGIC_out<-run_magic(data=t(RAW_DAT))
rownames(MAGIC_out)<-rownames(t(RAW_DAT))
colnames(MAGIC_out)<-colnames(t(RAW_DAT))
MAGIC_out<-t(MAGIC_out)
gr<-table(CONDITION)
M_magic<-SCnorm_runMAST(Data=MAGIC_out,NumCells = gr)

save(MAGIC_out,M_magic,file="E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/MAGIC_noDE_01_005_H1.RData")


####analysis######
load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")



#TrueMU<-rowData(SCE)$BaseGeneMean[-DROP]

Inputdat<-RAW_DAT

MedExp <- log(apply(Inputdat, 1, function(x) median(x[x != 0])))
MedExp[1:40]

# split into 6 equally sized groups:
grpnum <- 6
splitby <- sort(MedExp)
grps <- length(splitby)/grpnum
sreg <- split(splitby, ceiling(seq_along(splitby)/grps))

summary(sreg[[4]])

reee='MAST'



if(reee=='MAST'){

    #INPUT_LIST_temp<-list(apply(M_bayNorm_a10,1,median),M_scnorm$adjpval,M_scImpute$adjpval,M_RB$adjpval,apply(M_saver10,1,median),M_magic$adjpval,M_DCA$adjpval)
    INPUT_LIST_temp<-list(apply(M_bayNorm_a10,1,median),M_scnorm$adjpval,M_RB$adjpval,apply(M_saver10,1,median),M_magic$adjpval,M_DCA$adjpval)
}else if (reee=='TSTAT'){
    
    INPUT_LIST_temp<-list(T_mbayNorm,T_scnorm,T_scImpute,T_RB)
    
}

INPUT_LIST<-INPUT_LIST_temp
names(INPUT_LIST)<-names(INPUT_LIST_temp)

Gene_exp_gr<-seq(1,grpnum)

library(foreach)
thres<-0.05

norm_vec<-c('bayNorm','SCnorm','Scaling','SAVER','MAGIC','DCA')
BAR_MAST_l_list<-foreach(i=1:length(sreg))%do%{
    BAR_MAST_l<-lapply(INPUT_LIST,function(x){length(intersect(names(which(x<thres)),names(sreg[[i]])))})
    #BAR_MAST_l<-lapply(INPUT_LIST,function(x){length(intersect(names(which(x[,4]<thres)),names(sreg[[i]])))})
    
    return(BAR_MAST_l)
    
}

BAR_MAST_DAT<-foreach(i=1:length(BAR_MAST_l_list),.combine=rbind)%:%
    foreach(j=1:length(INPUT_LIST),.combine=rbind)%do%{
        temp<-BAR_MAST_l_list[[i]][[j]]
        temp<-c(temp,norm_vec[j],Gene_exp_gr[i])
        return(temp)
    }


BAR_MAST_DAT<-as.data.frame(BAR_MAST_DAT)
colnames(BAR_MAST_DAT)<-c('Number of detected DE genes','Normalization methods','Gene expression group')
BAR_MAST_DAT[,1]<-as.numeric(as.character(BAR_MAST_DAT[,1]))
BAR_MAST_DAT[,3]<-factor(BAR_MAST_DAT[,3],levels=unique(BAR_MAST_DAT[,3]))
BAR_MAST_DAT[,2]<-factor(BAR_MAST_DAT[,2],levels=unique(BAR_MAST_DAT[,2]))


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% norm_vec)]

library(ggplot2)
textsize<-14
DE_H1<-ggplot(data=BAR_MAST_DAT, aes(x=BAR_MAST_DAT[,3], y=BAR_MAST_DAT[,1], fill=BAR_MAST_DAT[,2])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=BAR_MAST_DAT[,1]), vjust=0, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "Gene expression group",y='Number of detected DE genes',fill='Normalization methods')+ggtitle("") +
    #scale_fill_brewer(palette="Paired")+
    scale_fill_manual(values=cbbPalette2)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
#dev.off()
DE_H1






######FC plot#####

DATA_list<-list(bayNorm=cbind(mbayNorm_out$Bay_mat_list$`Group 1`,mbayNorm_out$Bay_mat_list$`Group 2`),SCnorm=scnorm_out,Scaling=RB_out,SAVER=saver_out,MAGIC=MAGIC_out,DCA=DCA_out)

source("E:/RNAseqProject/Bacher__SCnorm_2016/FC_fun.R")
FC_out<-FC_fun(Inputdat=RAW_DAT,CONDITION=CONDITION,DATA_list,textsize=14,legend.key.size=1,colourval = cbbPalette2)
FC_out

library(gridExtra)
library(ggpubr)
library(cowplot)

qq<-plot_grid(DE_H1,FC_out + theme(legend.position="none"),ncol=1,nrow=2)

qq



# #####different DE method#####
# load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")
# source("E:/RNAseqProject/BAY_DEFUN.R")
# 
# baylist1<-list(mean=mbayNorm_out$Bay_mat_list$`Group 1`,mode=modebayNorm_out$Bay_mat_list$`Group 1`,array=abayNorm_out$Bay_array_list$`Group 1`)
# baylist2<-list(mean=mbayNorm_out$Bay_mat_list$`Group 2`,mode=modebayNorm_out$Bay_mat_list$`Group 2`,array=abayNorm_out$Bay_array_list$`Group 2`)
# 
# #DE_LIST<-BAY_DEFUN(baylist1,baylist2,CONDITION=CONDITION,NumCells=c(100,100))
# 
# load("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")
# source("E:/RNAseqProject/BAY_DEFUN.R")
# textsize<-12
# BAR_OUT_simH1<-BAR_DE_FUN(DE_list=DE_LIST,Inputdat=RAW_DAT,xlab='DE method_type of output from bayNorm')
# BAR_OUT_simH1
# 
# library(ggplot2)
# ggsave(filename="E:/RNAseqProject/Illustrator_bayNorm/SUP/SUP_SIMnoDE_H1.pdf",plot=BAR_OUT_simH1,width=8.2,height=5.5,units='in')



save.image("E:/RNAseqProject/SIMULATION/SIM_noDE_01_005_H1/SIM_noDE_01_005_H1.RData")
