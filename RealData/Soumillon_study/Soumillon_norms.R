load("E:/RNAseqProject/Soumillon_2014/Soumillon_2014.RData")


####bayNorm######
library(bayNorm)
library(scran)
sf<-computeSumFactors(as.matrix(D3_used3))
summary(sf/mean(sf)*0.03)
BETA<-sf/mean(sf)*0.03
names(BETA)<-colnames(D3_used3)

aBay_out<-bayNorm(Data=as.matrix(D3_used3),BETA_vec =BETA ,Conditions=CONDITION[colnames(D3_used3)],S=5,Prior_type = 'LL')



library(abind)
qq<-abind(aBay_out$Bay_array_list$`Group 1`,aBay_out$Bay_array_list$`Group 2`,along=2)
M_aBay<-SCnorm_runMAST3(Data=qq,NumCells = as.numeric(table(CONDITION[colnames(D3_used3)])))


####other norms############
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

#Scaling method####
RB<-t(t(D3_used3)/do.call(c,aBay_out$BETA))
M_RB<-SCnorm_runMAST3(Data=RB,NumCells = as.numeric(table(CONDITION[colnames(D3_used3)])))


#MAGIC####
library(Rmagic)
MAGIC_out<-run_magic(data=t(D3_used3))
MAGIC_out<-t(MAGIC_out)
rownames(MAGIC_out)<-rownames(D3_used3)
colnames(MAGIC_out)<-colnames(D3_used3)
M_MAGIC<-SCnorm_runMAST3(Data=MAGIC_out,NumCells = as.numeric(table(CONDITION[colnames(D3_used3)])))

#SCnorm#####
library(SCnorm)
SCnorm_out<-SCnorm(Data=D3_used3,Conditions=CONDITION[colnames(D3_used3)],ditherCounts = T)
M_SCnorm<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData,NumCells = as.numeric(table(CONDITION[colnames(D3_used3)])))



#scImpute#####
library(scImpute)
scimpute_fun(Data=D3_used3,Data_name='scImpute_out',labeled=T,labels=CONDITION[colnames(D3_used3)])
scImpute_out<-readRDS("E:/RNAseqProject/Soumillon_2014/scImpute_outscimpute_count.rds")


M_scimpute<-SCnorm_runMAST3(Data=scImpute_out,NumCells = as.numeric(table(CONDITION[colnames(D3_used3)])))


##SAVER####
library(SAVER)
saver_gr1<-saver(x=D3_used3[,CONDITION[colnames(D3_used3)]==1])
saver_gr2<-saver(x=D3_used3[,CONDITION[colnames(D3_used3)]==2])

saver_gr1_s<-sample.saver(saver_gr1,rep=5)
saver_gr2_s<-sample.saver(saver_gr2,rep=5)

library(abind)
saver_gr1_a<-abind(saver_gr1_s,along=3)
saver_gr2_a<-abind(saver_gr2_s,along=3)

sum(is.na(saver_gr1_a))
sum(is.na(saver_gr2_a))
saver_gr1_a[is.na(saver_gr1_a)]<-0
saver_gr2_a[is.na(saver_gr2_a)]<-0

qq<-abind(saver_gr1_a,saver_gr2_a,along=2)
dim(qq)
M_saver<-SCnorm_runMAST3(qq,NumCells=c(832,922))
M_SAVER<-SCnorm_runMAST3(Data=cbind(saver_gr1$estimate,saver_gr2$estimate),NumCells = c(832,922))




save.image(file="E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")