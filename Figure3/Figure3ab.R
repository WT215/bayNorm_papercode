#load("E:/RNAseqProject/PROJECT_scRNAseq/FIGURE_DROPOUT/DROPOUT_MAIN/Klein933_DIY.RData")
load('Klein_bayNorm.RData')
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

Real_data<-Real_Klein
length(which(rowMeans(Real_data)==0))
#sort the total counts of original Klein study data
qor<-sort(colSums(Real_data),decreasing = F)

#Pick the 100 cells with the smallest total counts and another 100 cells with the largest total counts to form two groups
Group1<-Real_data[,names(qor)[seq(1,100)]]
Group2<-Real_data[,rev(names(qor))[seq(1,100)]]
dim(Group1)

INPUT_DAT<-cbind(Group1,Group2)
length(which(rowSums(INPUT_DAT)==0))

INPUT_DAT<-INPUT_DAT[-which(rowSums(INPUT_DAT)==0),]


# boxplot(list(colSums(Group1),colSums(Group2)))
# library(bayNorm)
# BETA_out<-BetaFun(cbind(Group1,Group2),MeanBETA = 0.06)
# boxplot(list(BETA_out$BETA[colnames(Group1)],BETA_out$BETA[colnames(Group2)]))
# summary(BETA_out$BETA)

#use scran to estimate BETA
library(scran)
sf<-computeSumFactors(as.matrix(INPUT_DAT))
BETA_scran<-sf/mean(sf)*0.06
summary(BETA_scran)


#bayNorm####
library(bayNorm)
bayklein_de<-bayNorm(Data=INPUT_DAT,BETA_vec = BETA_scran,S=5,Conditions=c(rep(1,100),rep(2,100)),Prior_type = 'GG')

library(abind)
qq<-abind(bayklein_de$Bay_array_list,along=2)
M_bay<-SCnorm_runMAST3(Data=qq,NumCells = c(100,100))
length(which(apply(M_bay,1,median)<0.05))

mbayklein_de<-bayNorm_sup(Data=INPUT_DAT,BETA_vec = BETA_scran,S=1000,mean_version = T,Conditions=c(rep(1,100),rep(2,100)),PRIORS = bayklein_de$PRIORS_LIST)

save.image("E:/RNAseqProject/Klein2015/Klein_de.RData")
load("E:/RNAseqProject/Klein2015/Klein_de.RData")





####MAGIC#######
library(Rmagic)
magic_out<-run_magic(data=t(INPUT_DAT))
magic_out<-t(magic_out)
colnames(magic_out)<-colnames(INPUT_DAT)
M_magic<-SCnorm_runMAST3(Data=magic_out,NumCells = c(100,100))
save(magic_out,M_magic,file="E:/RNAseqProject/Klein2015/MAGIC_out.RData")
length(which(M_magic$adjpval<0.05))

#####RB########
RB_out<-t(t(INPUT_DAT/do.call(c,bayklein_de$BETA)))
M_RB<-SCnorm_runMAST3(Data=RB_out,NumCells = c(100,100))
save(RB_out,M_RB,file="E:/RNAseqProject/Klein2015/RB_out.RData")

###scImpute#####
library(scImpute)
scimpute_fun(Data=INPUT_DAT,Data_name='scImpute_klein',Kcluster=1)
scImpute_out<-readRDS("E:/RNAseqProject/Klein2015/scImpute_kleinscimpute_count.rds")
M_scImpute<-SCnorm_runMAST3(Data=scImpute_out,NumCells = c(100,100))
save(scImpute_out,M_scImpute,file="E:/RNAseqProject/Klein2015/scImpute_out.RData")

length(which(M_scImpute$adjpval<0.05))
###scnorm####
library(SCnorm)
SCnorm_out<-SCnorm(Data=INPUT_DAT,Conditions = c(rep(1,100),rep(2,100)),ditherCounts = T)
M_SCnorm<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData,NumCells = c(100,100))

save(SCnorm_out,M_SCnorm,file="E:/RNAseqProject/Klein2015/SCnorm_out.RData")

#DCA#####
load("E:/RNAseqProject/Klein2015/DCA_klein.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_klein,NumCells = c(100,100))
length(which(M_DCA$adjpval<0.05))
save(DCA_klein,M_DCA,file="E:/RNAseqProject/Klein2015/DCA_klein.RData")


####SAVER##########
library(SAVER)
saver_out<-saver(x=INPUT_DAT)
samplesave<-sample.saver(x=saver_out,rep=5)
library(abind)
saver_array<-abind(samplesave,along=3)
dim(saver_array)
saver_array[is.na(saver_array)]<-0
M_saver<-SCnorm_runMAST3(Data=saver_array, NumCells=c(100,100))
save(saver_out,saver_array,M_saver,file="E:/RNAseqProject/Klein2015/saver_out.RData")

length(which(apply(M_saver,1,median)<0.05))


#######begin analysis######
load("E:/RNAseqProject/Klein2015/Klein_de.RData")



load("E:/RNAseqProject/Klein2015/MAGIC_out.RData")
load("E:/RNAseqProject/Klein2015/RB_out.RData")
load("E:/RNAseqProject/Klein2015/scImpute_out.RData")
load("E:/RNAseqProject/Klein2015/SCnorm_out.RData")
load("E:/RNAseqProject/Klein2015/DCA_klein.RData")
load("E:/RNAseqProject/Klein2015/saver_out.RData")


DATA_list<-list(bayNorm=do.call(cbind,mbayklein_de$Bay_mat_list),SCnorm=SCnorm_out@metadata$NormalizedData,Scaling=RB_out,SAVER=saver_out$estimate,MAGIC=magic_out,DCA=DCA_klein)

norm_vec<-c('bayNorm','SCnorm','Scaling','SAVER','MAGIC','DCA')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% norm_vec)][norm_vec]

source("E:/RNAseqProject/Bacher__SCnorm_2016/FC_fun.R")
FC_H1<-FC_fun(Inputdat=INPUT_DAT,CONDITION=c(rep(1,100),rep(2,100)),DATA_list,textsize=10,legend.key.size=0.5,colourval=cbbPalette2)
FC_H1






# MAST results############
Inputdat<-INPUT_DAT

MedExp <- log(apply(Inputdat, 1, function(x) median(x[x != 0])))
summary(apply(Inputdat, 1, function(x) median(x[x != 0])))
length(which(rowSums(Inputdat)==0))


hist(exp(MedExp[names(sreg$`6`)]),breaks=200)
summary(exp(MedExp[names(sreg$`6`)]))


# split into 4 equally sized groups:
grpnum <- 6
splitby <- sort(MedExp)
grps <- length(splitby)/grpnum
sreg <- split(splitby, ceiling(seq_along(splitby)/grps))
lapply(sreg,length)


CONDITION<-c(rep(1,100),rep(2,100))

INPUT_LIST_temp<-list(R_MAST_bay_adj=apply(M_bay,1,median),R_MAST_scnorm_adj=M_SCnorm$adjpval,Scaling=M_RB$adjpval,SAVER=apply(M_saver,1,median),MAGIC=M_magic$adjpval,DCA=M_DCA$adjpval)




INPUT_LIST<-INPUT_LIST_temp
names(INPUT_LIST)<-names(INPUT_LIST_temp)

Gene_exp_gr<-seq(1,grpnum)

library(foreach)
thres<-0.05

norm_vec<-c('bayNorm','SCnorm','Scaling','SAVER','MAGIC','DCA')
BAR_MAST_l_list<-foreach(i=1:length(sreg))%do%{
    BAR_MAST_l<-lapply(INPUT_LIST,function(x){length(intersect(names(which(x<thres)),names(sreg[[i]])))})
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
textsize<-10
useee<-BAR_MAST_DAT$`Normalization methods`!='scImpute'
DE_H1<-ggplot(data=BAR_MAST_DAT[useee,], aes(x=BAR_MAST_DAT[useee,3], y=BAR_MAST_DAT[useee,1], fill=BAR_MAST_DAT[useee,2])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=BAR_MAST_DAT[useee,1]), vjust=0, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "Gene expression group",y='Number of detected DE genes',fill='Normalization methods')+ggtitle("Subset of Klein study") +
    #scale_fill_brewer(palette="Paired")+
    scale_fill_manual(values=cbbPalette2)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
#dev.off()
DE_H1



####Fig 3(a)-(b)#######
source(file="E:/RNAseqProject/MANY_SAVE_PATH.R")
library(gridExtra)
library(ggpubr)
library(cowplot)
#qq<-ggarrange(DE_H1 ,MSE_H1,ncol=2,nrow=1,common.legend = TRUE, legend="top")
qq<-plot_grid(DE_H1 + theme(legend.position="none"),FC_H1,ncol=2,nrow=1)
#draw_grob(get_legend(DE_H1), 0.55, 0.75, 1/3, 0.5)
qq
ggsave(FIGURE_2_PATH_fun("/FIG2_klein_5samplestr_scranbeta.pdf"),plot=qq,width =11, height =5,units='in')
