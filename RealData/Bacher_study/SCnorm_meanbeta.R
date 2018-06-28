###H1#####
load("RAW_INITIATE.RData")


H1_data_comb<-as.matrix(cbind(H1_p24,H1_p96)[whichg_H1,])

colSums(H1_data_comb)


CONDITION_H1<-c(rep(1,dim(H1_p24)[2]),rep(2,dim(H1_p96)[2]))

#5% mean BETA######
ERCC_BETA<-colSums(cbind(ERCC_H1_p24/20,ERCC_H1_p96/10))
ERCC_BETA_005<-ERCC_BETA/mean(ERCC_BETA)*0.05
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")


library(bayNorm)
mBAY_H1_005<-bayNorm(Data=H1_data_comb,BETA_vec=ERCC_BETA_005,S=1000,parallel=T,NCores=5,FIX_MU = T,GR=F,Conditions=CONDITION_H1,BB_SIZE = T,mode_version = F,mean_version = T,UMI_sffl=c(20,10),Prior_type = 'GG',verbose = T)

M_mBAY_H1_005<-SCnorm_runMAST(Data=cbind(mBAY_H1_005$Bay_mat_list$`Group 1`,mBAY_H1_005$Bay_mat_list$`Group 2`),NumCells=c(92,92))



save(H1_data_comb,mBAY_H1_005,M_mBAY_H1_005,T_mBAY_H1_005,file="E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/meanBay/H1_BAY_esf_005.RData")


####20% mean BETA#####
ERCC_BETA_02<-ERCC_BETA/mean(ERCC_BETA)*0.2

source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")


library(bayNorm)
mBAY_H1_02<-bayNorm(Data=H1_data_comb,BETA_vec=ERCC_BETA_02,S=1000,parallel=T,NCores=5,FIX_MU = T,GR=F,Conditions=CONDITION_H1,BB_SIZE = T,mode_version = F,mean_version = T,UMI_sffl=c(20,10),Prior_type = 'GG',verbose = T)

M_mBAY_H1_02<-SCnorm_runMAST(Data=cbind(mBAY_H1_02$Bay_mat_list$`Group 1`,mBAY_H1_02$Bay_mat_list$`Group 2`),NumCells=c(92,92))



save(H1_data_comb,mBAY_H1_02,M_mBAY_H1_02,T_mBAY_H1_02,file="E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/meanBay/H1_BAY_esf_02.RData")





#######analysis#####
load("H1_many_normalizations.RData")

load("E:/RNAseqProject/Bacher__SCnorm_2016/RAW_INITIATE.RData")
load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/meanBay/H1_BAY_esf_02.RData")
load("E:/RNAseqProject/Bacher__SCnorm_2016/scaled_BAY/meanBay/H1_BAY_esf_005.RData")




H1_data_comb<-cbind(H1_p24,H1_p96)[whichg_H1,]

# Conditional median over all cells.
Inputdat<-H1_data_comb

#The following code was kindly provided by Bacher (The author of SCnorm)
MedExp <- log(apply(Inputdat, 1, function(x) median(x[x != 0])))
# split into 4 equally sized groups:
grpnum <- 6
splitby <- sort(MedExp)
grps <- length(splitby)/grpnum
sreg <- split(splitby, ceiling(seq_along(splitby)/grps))

CONDITION=CONDITION_H1



reee='MAST'
INPUT_LIST_temp<-list(M_mBAY_H1$adjpval,M_mBAY_H1_005$adjpval,M_mBAY_H1_02$adjpval)

INPUT_LIST<-INPUT_LIST_temp

for(i in 1:length(INPUT_LIST)){
    names(INPUT_LIST[[i]])<-names(M_mBAY_H1$adjpval)
}

names(INPUT_LIST)<-names(INPUT_LIST_temp)

Gene_exp_gr<-seq(1,grpnum)

library(foreach)
#Set the threshold for the adjusted P-values to be 0.05
thres<-0.05

norm_vec<-c('bayNorm_01','bayNorm_005','bayNorm_02')
BAR_MAST_l_list<-foreach(i=1:length(sreg))%do%{
    BAR_MAST_l<-lapply(INPUT_LIST,function(x){length(intersect(names(which(x<thres)),names(sreg[[i]])))})
    return(BAR_MAST_l)
    
}

BAR_MAST_names_list<-foreach(i=1:length(sreg))%do%{
    BAR_MAST_l<-lapply(INPUT_LIST,function(x){intersect(names(which(x<thres)),names(sreg[[i]]))})
    return(BAR_MAST_l)
    
}
for(i in 1:4){
    names(BAR_MAST_names_list[[i]])<-c('bayNorm_01','bayNorm_005','bayNorm_02')
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

library(ggplot2)
textsize<-6
DE_H1<-ggplot(data=BAR_MAST_DAT, aes(x=BAR_MAST_DAT[,3], y=BAR_MAST_DAT[,1], fill=BAR_MAST_DAT[,2])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_text(aes(label=BAR_MAST_DAT[,1]), vjust=1.6, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "Gene expression group",y='Number of detected DE genes',fill='Normalization methods')+ggtitle("") +
    scale_fill_brewer(palette="Paired")+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='top')
#dev.off()
DE_H1


######FC plot: Fig S23 (a)#####
DATA_list<-list(bayNorm_01=cbind(mBAY_H1$Bay_mat_list$`Group 1`,mBAY_H1$Bay_mat_list$`Group 2`),bayNorm_005=do.call(cbind,mBAY_H1_005$Bay_mat_list),bayNorm_02=do.call(cbind,mBAY_H1_02$Bay_mat_list))
names(DATA_list)<-c("bayNorm_10%" , "bayNorm_5%", "bayNorm_20%" )

source("E:/RNAseqProject/Bacher__SCnorm_2016/FC_fun.R")
FC_H1<-FC_fun(Inputdat=cbind(H1_p24,H1_p96)[whichg_H1,],CONDITION=CONDITION_H1,DATA_list,textsize=14,legend.key.size=1,colourval=c(1,2,3))
FC_H1

