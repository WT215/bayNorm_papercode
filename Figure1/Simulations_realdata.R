#load necessary functions########
source("E:/RNAseqProject/QQsim_v2_SingleCellExperiment.R")
sm <- function(x, y, x.log = FALSE,n.bins = 25){
    if(x.log){ 
        brks <- unique(quantile(x, probs = seq(0,1,len=25))) 
    } else {
        brks <- 2^unique(quantile(log2(x), probs = seq(0,1,len=n.bins))) 
    }
    mids <- (brks[-1] + brks[-length(brks)] )/ 2
    x.in <- cut(x, breaks = brks, include.lowest = TRUE)
    m <- tapply(y, x.in, mean)
    fit = lm(y~x)
    l <- predict(fit, newdata = data.frame(x  = mids))
    dat <- data.frame(x=mids, y = m, n = as.numeric(table(x.in)), pred = l)
}


dropoutfun<-function(data){
    xx<-rowMeans(data)
    yy<-apply(data,1,function(x){length(which(x==0))/length(x)})
    qq<-cbind(xx,yy)
    return(qq)
}
SIM_FUN<-function(DATA,MU,SIZE,BETA)
{
    
    nCells<-dim(DATA)[2]
    nGenes<-dim(DATA)[1]
    
    GeneMean_mat<-matrix(MU,ncol=nCells ,nrow=nGenes,byrow=F)
    
    one_bcv2<-SIZE
    
    Gamma_Means_mat <- matrix(rgamma(nGenes * nCells, shape =one_bcv2, scale = GeneMean_mat * (1/one_bcv2)),nrow = nGenes, ncol = nCells)
    
    true.counts <- matrix(rpois(nGenes * nCells, lambda =  Gamma_Means_mat ),nrow = nGenes, ncol = nCells)
    rownames(true.counts)<-rownames(DATA)
    colnames(true.counts)<-colnames(DATA)
    
    downsample.counts <-bayNorm::DownSampling(true.counts,BETA)
    
    rownames(downsample.counts)<-rownames(DATA)
    colnames(downsample.counts)<-colnames(DATA)
    
    return(list(true.counts=true.counts,downsample.counts=downsample.counts))
}


#######REAL DATA 1: Klein study####
# #Import Klein mouse cells 933 cells
GSM1599494_ES_d0_main <- read.csv("E:/RNAseqProject/NEWPROJECT_PAPERS/Droplet Barcoding for Single-Cell Transcriptomics/allcsv/GSM1599494_ES_d0_main.csv",header=FALSE)
Whole_dat<-as.data.frame(GSM1599494_ES_d0_main)
rownames(Whole_dat)<-Whole_dat[,1]
Whole_dat<-Whole_dat[,-1]
rm(GSM1599494_ES_d0_main)

#estimate BETA for bayNorm
sf<-apply(Whole_dat,2,mean,trim=0.01)
BETA_Klein<-sf/mean(sf)*0.06
library(bayNorm)
Real_Klein<-Whole_dat
rm(Whole_dat)

#begin bayNorm
bay_Klein<-bayNorm(Data=Real_Klein,BETA_vec = BETA_Klein,mean_version = T)



#run "Binomial_bayNorm" simulation protocol
BaySim_Kelin<-SIM_FUN(DATA=Real_Klein,MU=bay_Klein$PRIORS$MME_prior$MME_MU,SIZE=bay_Klein$PRIORS$MME_SIZE_adjust,BETA=bay_Klein$BETA)

#run "Binomial" simulation protocol
Sim_List_Input<-trytt2(counts=Real_Klein,inputBeta = NULL,inputMeanBeta = 0.06,inputtrim = 0.01)

#run Splatter
library(splatter)
splatter_klein_params <- splatEstimate(as.matrix(Real_data))
splatter_klein_sim<-splatSimulate(splatter_klein_params)

#Prepare SCE lists, then we will put it into a function for producing comparison plots
SCElist_Klein2<-list(
    SCElist_Klein$Real,
    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_Kelin$downsample.counts)),
    Binomial=SingleCellExperiment(assays=list(counts=as.matrix(Sim_List_Input$SCE@assays@.xData$data$counts))),
    Splatter=splatter_klein_sim)



#####main Fig 1: Klein####### 
point.size<-1
point.alpha<-0.4
linewidth<-1
linewidth.exp<-1.5

textsize<-10
legendpointsize=4
legend_key_size=0.8



source("DROPOUT_FUN.r")
meanvar_klein<-plot_MEANVAR_v2(sces=SCElist_Klein2,MAIN='Klein et al study',CAPTION='',legend.position=c(0.8,0.8))

library(gridExtra)
library(ggpubr)

qq<-ggarrange(
    
    as_ggplot(meanvar_klein$mean.var),
    
    as_ggplot(meanvar_klein$mean.zeros),
    
    
    as_ggplot( meanvar_klein$z.gene),
    
    as_ggplot(meanvar_klein$z.cell),
    
    ncol=2,nrow=2,common.legend = T, legend="top")
qq






#Sup figure for Klein#####
point.size<-1
point.alpha<-0.4
linewidth<-0.5
linewidth.exp<-1

textsize<-6
legendpointsize=4
legend_key_size=0.8
legend.position='top'
hline_alpha<-0.5
hline_size<-0.8


source("DROPOUT_FUN.r")
meanvar_klein<-plot_MEANVAR(sces=SCElist_Klein2,MAIN='Klein et al study',CAPTION='',legend.position=c(0.8,0.85))
meanvar_klein_diff<-plot_MEANVAR_diff(sces=SCElist_Klein2,MAIN='',CAPTION='')


library(gridExtra)
library(ggpubr)

qq<-ggarrange(
    
    meanvar_klein$mean.var,
    meanvar_klein_diff$mean.var,
    
    
    meanvar_klein$mean.zeros,
    meanvar_klein_diff$mean.zeros,
    
    
    meanvar_klein$z.gene,
    meanvar_klein_diff$z.gene,
    
    meanvar_klein$z.cell,
    meanvar_klein_diff$z.cell,
    
    ncol=2,nrow=4,common.legend = F)
qq





#######REAL DATA 2: Tung study (individual NA19098)####
rm(list=ls())


load("E:/RNAseqProject/Illustrator_bayNorm/bayNorm_papercode/RealData/Tung_study/Load_Tung.RData")
load('Tung_norms.RData')

library(SingleCellExperiment)
library(bayNorm)

#run "Binomial_bayNorm" simulation protocol
BaySim_TungN1<-SIM_FUN(DATA=N1_DAT,MU=bayNorm_mean_N1$PRIORS$MME_prior$MME_MU,SIZE=bayNorm_mean_N1$PRIORS$MME_SIZE_adjust,BETA=bayNorm_mean_N1$BETA)


#run "Binomial" simulation protocol
TungN1_Binom<-trytt2(N1_DAT,inputMeanBeta=1,inputtrim=0.01,inputBeta=efficiency[colnames(N1_DAT)])

#run Splatter
library(splatter)
TungN1_est_splatter<-splatEstimate(counts=N1_DAT)
TungN1_sim_splatter<-splatSimulate(params=TungN1_est_splatter)

#Prepare SCE lists, then we will put it into a function for producing comparison plots
SCElist_Tung2<-list(
    SCElist_Tung$Real,
    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_TungN1$downsample.counts)),
    Binomial=SingleCellExperiment(assays=list(counts=TungN1_Binom$SCE@assays@.xData$data$counts)),
    Splatter=SingleCellExperiment(assays=list(counts=TungN1_sim_splatter@assays@.xData$data$counts)))
names(SCElist_Tung2)[1]<-'Real data'



names(SCElist_Tung2)<-c('Real data','Binomial_bayNorm','Binomial','Splatter')


point.size<-1
point.alpha<-0.4
linewidth<-0.5
linewidth.exp<-1
textsize<-6
legendpointsize=4
legend_key_size=0.8
legend.position='top'
hline_alpha<-0.5
hline_size<-0.8
library(SingleCellExperiment)
source("DROPOUT_FUN.r")
meanvar_Tung<-plot_MEANVAR(sces=SCElist_Tung2,MAIN='Data from Tung case study (NA19098)',CAPTION='',legend.position = c(0.8,0.85))
meanvar_Tung_diff<-plot_MEANVAR_diff(sces=SCElist_Tung2,MAIN='',CAPTION='')




library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    meanvar_Tung$mean.var,
    meanvar_Tung_diff$mean.var,
    
    meanvar_Tung$mean.zeros,
    meanvar_Tung_diff$mean.zeros,
    
    
    
    meanvar_Tung$z.gene,
    meanvar_Tung_diff$z.gene,
    
    meanvar_Tung$z.cell,
    meanvar_Tung_diff$z.cell,
    
    ncol=2,nrow=4,common.legend =F)

qq


save.image('Tung_N1_SIM')



#######REAL DATA 3: Tung study (individual NA19101)####
load("E:/RNAseqProject/Illustrator_bayNorm/bayNorm_papercode/RealData/Tung_study/Load_Tung.RData")
load('Tung_norms.RData')

library(SingleCellExperiment)
library(bayNorm)

#run "Binomial_bayNorm" simulation protocol
BaySim_TungN2<-SIM_FUN(DATA=N2_DAT,MU=bayNorm_mean_N2$PRIORS$MME_prior$MME_MU,SIZE=bayNorm_mean_N2$PRIORS$MME_SIZE_adjust,BETA=bayNorm_mean_N2$BETA)


#run "Binomial" simulation protocol
TungN2_Binom<-trytt2(N2_DAT,inputMeanBeta=1,inputtrim=0.01,inputBeta=bayNorm_mean_N2$BETA)

#run Splatter
library(splatter)
TungN2_est_splatter<-splatEstimate(counts=N2_DAT)
TungN2_sim_splatter<-splatSimulate(params=TungN2_est_splatter)

#Prepare SCE lists, then we will put it into a function for producing comparison plots
SCElist_TungN2<-list(
    SingleCellExperiment(assays=list(counts=N2_DAT)),
    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_TungN2$downsample.counts)),
    Binomial=SingleCellExperiment(assays=list(counts=TungN2_Binom$SCE@assays@.xData$data$counts)),
    Splatter=SingleCellExperiment(assays=list(counts=TungN2_sim_splatter@assays@.xData$data$counts)))
names(SCElist_TungN2)[1]<-'Real data'

source("DROPOUT_FUN.r")
meanvar_TungN2<-plot_MEANVAR(sces=SCElist_TungN2,MAIN='Data from Tung case study (NA19101)',CAPTION='',legend.position = c(0.8,0.85))
meanvar_TungN2_diff<-plot_MEANVAR_diff(sces=SCElist_TungN2,MAIN='',CAPTION='')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    meanvar_TungN2$mean.var,
    meanvar_TungN2_diff$mean.var,
    
    meanvar_TungN2$mean.zeros,
    meanvar_TungN2_diff$mean.zeros,
    
    meanvar_TungN2$z.gene,
    meanvar_TungN2_diff$z.gene,
    
    meanvar_TungN2$z.cell,
    meanvar_TungN2_diff$z.cell,
    
    ncol=2,nrow=4,common.legend = F)

qq


#######REAL DATA 4: Tung study (individual NA19239)####
load("E:/RNAseqProject/Illustrator_bayNorm/bayNorm_papercode/RealData/Tung_study/Load_Tung.RData")
load('Tung_norms.RData')

library(SingleCellExperiment)
library(bayNorm)

#run "Binomial_bayNorm" simulation protocol
BaySim_TungN3<-SIM_FUN(DATA=N3_DAT,MU=bayNorm_mean_N3$PRIORS$MME_prior$MME_MU,SIZE=bayNorm_mean_N3$PRIORS$MME_SIZE_adjust,BETA=bayNorm_mean_N3$BETA)


#run "Binomial" simulation protocol
TungN3_Binom<-trytt2(N3_DAT,inputMeanBeta=1,inputtrim=0.01,inputBeta=bayNorm_mean_N3$BETA)


#run Splatter
library(splatter)
TungN3_est_splatter<-splatEstimate(counts=N3_DAT)
TungN3_sim_splatter<-splatSimulate(params=TungN3_est_splatter)

#Prepare SCE lists, then we will put it into a function for producing comparison plots
SCElist_TungN3<-list(
    SingleCellExperiment(assays=list(counts=N3_DAT)),
    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_TungN3$downsample.counts)),
    Binomial=SingleCellExperiment(assays=list(counts=TungN3_Binom$SCE@assays@.xData$data$counts)),
    Splatter=SingleCellExperiment(assays=list(counts=TungN3_sim_splatter@assays@.xData$data$counts)))
names(SCElist_TungN3)[1]<-'Real data'


point.size<-1
point.alpha<-0.4
linewidth<-0.5
linewidth.exp<-1

textsize<-6
legendpointsize=4
legend_key_size=0.8
legend.position='top'
hline_alpha<-0.5
hline_size<-0.8
source("DROPOUT_FUN.r")
meanvar_TungN3<-plot_MEANVAR(sces=SCElist_TungN3,MAIN='Data from Tung case study (NA19239)',CAPTION='',legend.position = c(0.8,0.85))
meanvar_TungN3_diff<-plot_MEANVAR_diff(sces=SCElist_TungN3,MAIN='',CAPTION='')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    meanvar_TungN3$mean.var,
    meanvar_TungN3_diff$mean.var,
    
    meanvar_TungN3$mean.zeros,
    meanvar_TungN3_diff$mean.zeros,
    
    meanvar_TungN3$z.gene,
    meanvar_TungN3_diff$z.gene,
    
    meanvar_TungN3$z.cell,
    meanvar_TungN3_diff$z.cell,
    
    ncol=2,nrow=4,common.legend = F)

qq


#######REAL DATA 5: Torre study####
rm(list=ls())

load("Torre_many_normalizations.RData")


#run "Binomial_bayNorm" simulation protocol
BaySim_Torre<-SIM_FUN(DATA=Torre_drop_sub,MU=bay_out$PRIORS$MME_prior$MME_MU,SIZE=bay_out$PRIORS$MME_SIZE_adjust,BETA=bay_out$BETA)

#run "Binomial" simulation protocol
Torre_Binom<-trytt2(Torre_drop_sub,inputMeanBeta=1,inputtrim=0.01,inputBeta=bay_out$BETA)

#run Splatter
library(splatter)
Torre_est_splatter<-splatEstimate(counts=Torre_drop_sub)
Torre_sim_splatter<-splatSimulate(params=Torre_est_splatter)




SCElist_Torre<-list(Raw=SingleCellExperiment(assays=list(counts=Torre_drop_sub)),
                    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_Torre$downsample.counts)),
                    Binomial=SingleCellExperiment(assays=list(counts=Torre_Binom$SCE@assays@.xData$data$counts)),
                    Splatter=SingleCellExperiment(assays=list(counts=Torre_sim_splatter@assays@.xData$data$counts)))
names(SCElist_Torre)[1]<-'Real data'




source("DROPOUT_FUN.r")
meanvar_Torre<-plot_MEANVAR(sces=SCElist_Torre,MAIN='Data from Torre case study',CAPTION='',legend.position=c(0.8,1.10))
meanvar_Torre_diff<-plot_MEANVAR_diff(sces=SCElist_Torre,MAIN='',CAPTION='')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    meanvar_Torre$mean.var,
    meanvar_Torre_diff$mean.var,
    
    meanvar_Torre$mean.zeros,
    meanvar_Torre_diff$mean.zeros,
    
    meanvar_Torre$z.gene,
    meanvar_Torre_diff$z.gene,
    
    meanvar_Torre$z.cell,
    meanvar_Torre_diff$z.cell,
    
    ncol=2,nrow=4,common.legend = F)

qq


save.image('Torre_sim.RData')

#######REAL DATA 6: Bacher study (H1_P24 cells)####
load("H1_many_normalizations.RData")
library(splatter)
library(SingleCellExperiment)
library(bayNorm)
source("E:/RNAseqProject/QQsim_v2_SingleCellExperiment.R")


geneused<-rownames(mBAY_H1$Bay_mat_list$`Group 1`)
#run "Binomial_bayNorm" simulation protocol
Binomial_bayNorm<-SIM_FUN(DATA=H1_p24[geneused,],MU=mBAY_H1$PRIORS_LIST$`Group 1`$MME_prior$MME_MU,SIZE=mBAY_H1$PRIORS_LIST$`Group 1`$MME_SIZE_adjust,BETA=mBAY_H1$BETA$`Group 1`)



#run Splatter
count_temp<-round(as.matrix(H1_p24[geneused,]/20))
library(splatter)
spla_est<-splatEstimate(counts=count_temp)
spla_sim<-splatSimulate(spla_est)

#run "Binomial" simulation protocol
system.time(H1p24_bay_sim<-trytt2(counts=count_temp,inputMeanBeta=1,inputtrim=0.01,inputBeta=mBAY_H1$BETA$`Group 1`[colnames(count_temp)]))

sce_list<-list(real=SingleCellExperiment(assays=list(counts=count_temp)),Binomial_bayNorm=SingleCellExperiment(assays=list(counts=Binomial_bayNorm$downsample.counts)),Binomial=SingleCellExperiment(assays=list(counts=H1p24_bay_sim$SCE@assays$data$counts)),Splatter=spla_sim)
names(sce_list)[1]<-'Scaled and rounded real data'

source("DROPOUT_FUN.r")
point.size<-1
point.alpha<-0.4
linewidth<-0.5
linewidth.exp<-1

textsize<-6
legendpointsize=4
legend_key_size=0.8
legend.position='top'
hline_alpha<-0.5
hline_size<-0.8
library(ggplot2)
re1<-plot_MEANVAR(sce_list,MAIN='Data from Bacher study (H1_P24 cells)',CAPTION='',legend.position = c(0.7,0.8))
re2<-plot_MEANVAR_diff(sce_list,MAIN='',CAPTION='')


library(gridExtra)
library(ggpubr)

qq<-ggarrange(
    
    re1$mean.var,
    re2$mean.var,
    
    re1$mean.zeros,
    re2$mean.zeros,
    
    
    
    re1$z.gene,
    re2$z.gene,
    
    re1$z.cell,
    re2$z.cell,
    
    ncol=2,nrow=4,common.legend = F)
qq










######prepare for MA plot####
library(ggplot2)
library(SummarizedExperiment)
library(SingleCellExperiment)
source("DROPOUT_FUN.r")

load("Torre_many_normalizations.RData")
bay_out_torre<-bay_out


#SCnorm
load("H1_many_normalizations.RData")

geneused<-rownames(mBAY_H1$Bay_mat_list$`Group 1`)
Binomial_bayNorm_H1P24<-SIM_FUN(DATA=H1_p24[geneused,],MU=mBAY_H1$PRIORS_LIST$`Group 1`$MME_prior$MME_MU,SIZE=mBAY_H1$PRIORS_LIST$`Group 1`$MME_SIZE_adjust,BETA=mBAY_H1$BETA$`Group 1`)

#Patel
load('Patel2014_bay_out.RData')
bay_out_patel<-bay_out
MU<-bay_out_patel$PRIORS$MME_prior$MME_MU
SIZE<-bay_out_patel$PRIORS$MME_SIZE_adjust
Binomial_bayNorm_patel<-SIM_FUN(DATA=Inputdat,MU=MU,SIZE=SIZE,BETA=bay_out_patel$BETA)



##begin MA plot####
par(mfrow=c(3,2))
###MA_Torre
SIM_TrueMU<-bay_out_torre$PRIORS$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(BaySim_Torre$downsample.counts)/bay_out_torre$BETA))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Torre study (UMI)',col=8,ylim=c(-0.2,0.2))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
###MA_TUNG N1
SIM_TrueMU<-bayNorm_mean_N1$PRIORS$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(BaySim_TungN1$downsample.counts)/bayNorm_mean_N1$BETA))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Tung study (UMI, NA19098)',col=8,ylim=c(-0.2,0.2))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)


###MA_TUNG N2
SIM_TrueMU<-bayNorm_mean_N2$PRIORS$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(BaySim_TungN2$downsample.counts)/bayNorm_mean_N2$BETA))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Tung study (UMI, NA19101)',col=8,ylim=c(-0.2,0.2))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)

###MA_TUNG N3
SIM_TrueMU<-bayNorm_mean_N3$PRIORS$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(BaySim_TungN3$downsample.counts)/bayNorm_mean_N3$BETA))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Tung study (UMI, NA19239)',col=8,ylim=c(-0.2,0.2))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)



#Bacher
SIM_TrueMU<-mBAY_H1$PRIORS_LIST$`Group 1`$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(Binomial_bayNorm_H1P24$downsample.counts)/mBAY_H1$BETA$`Group 1`))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Bacher study (non UMI, H1_P24 cells)',col=8,ylim=c(-0.5,0.5))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)

#Patel
SIM_TrueMU<-bay_out_patel$PRIORS$MME_prior$MME_MU
SIM_ESTMU<-rowMeans(t(t(Binomial_bayNorm_patel$downsample.counts)/bay_out_patel$BETA))
yy<-log2(SIM_ESTMU+1)-log2(SIM_TrueMU+1)
xx<-(log2(SIM_ESTMU+1)+log2(SIM_TrueMU+1))/2
plot(xx,yy,xlab='A=(log2(EST MU)+log2(TRUE MU))/2',ylab='M=log2(EST MU)-log2(TRUE MU)',pch='.',main='Patel study (non UMI)',col=8,ylim=c(-0.5,0.5))
abline(h = 0, col = 2, lwd = 1)
smDat <- sm(x = xx, y = yy,n.bins =50)
lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)


