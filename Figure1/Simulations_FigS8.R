######################### load functions##########
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


SIM_FUN_splatter<-function(DATA,MU,SIZE,BETA,bcv.df,bcv.common)
{
    
    nCells<-dim(DATA)[2]
    nGenes<-dim(DATA)[1]
    
    GeneMean_mat<-matrix(MU,ncol=nCells ,nrow=nGenes,byrow=F)
    
    
    GeneMean_mat2<-t(t(GeneMean_mat)*BETA)
    
    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(GeneMean_mat2))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    }else{
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(GeneMean_mat2)))
    }
    
    
    one_bcv2<-1/(bcv^2)
    one_bcv2[one_bcv2==0]=min(one_bcv2[one_bcv2>0])
    #one_bcv2<-SIZE
    
    
    
    Gamma_Means_mat <- matrix(rgamma(nGenes * nCells, shape =one_bcv2, scale = GeneMean_mat2 * (1/one_bcv2)),nrow = nGenes, ncol = nCells)
    
    true.counts <- matrix(rpois(nGenes * nCells, lambda =  Gamma_Means_mat ),nrow = nGenes, ncol = nCells)
    rownames(true.counts)<-rownames(DATA)
    colnames(true.counts)<-colnames(DATA)
    
    #downsample.counts <-bayNorm::DownSampling(true.counts,BETA)
    
    #rownames(downsample.counts)<-rownames(DATA)
    #colnames(downsample.counts)<-colnames(DATA)
    
    return(list(true.counts=true.counts))
}



SIM_FUN_saver<-function(DATA,saverdat,BETA)
{
    
    nCells<-dim(DATA)[2]
    nGenes<-dim(DATA)[1]
    
    #GeneMean_mat<-matrix(MU,ncol=nCells ,nrow=nGenes,byrow=F)
    
    
    Gamma_Means_mat<-t(t(saverdat)*BETA)
    
    
    #one_bcv2<-SIZE
    
    #Gamma_Means_mat <- matrix(rgamma(nGenes * nCells, shape =one_bcv2, scale = GeneMean_mat2 * (1/one_bcv2)),nrow = nGenes, ncol = nCells)
    
    true.counts <- matrix(rpois(nGenes * nCells, lambda =  Gamma_Means_mat ),nrow = nGenes, ncol = nCells)
    rownames(true.counts)<-rownames(DATA)
    colnames(true.counts)<-colnames(DATA)
    
    #downsample.counts <-bayNorm::DownSampling(true.counts,BETA)
    
    #rownames(downsample.counts)<-rownames(DATA)
    #colnames(downsample.counts)<-colnames(DATA)
    
    return(list(true.counts=true.counts))
}
####Klein####

load('Klein_bayNorm.RData')

library(ggplot2)
source("DROPOUT_FUN.r")


library(SingleCellExperiment)
real=SingleCellExperiment(assays=list(counts=as.matrix(Real_Klein)))




trim001_012<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta = 0.12,inputBeta=NULL,inputtrim = 0.01)
trim001_006<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta = 0.06,inputBeta=NULL,inputtrim = 0.01)
trim001_003<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta = 0.03,inputBeta=NULL,inputtrim = 0.01)

library(scran)
scransf<-computeSumFactors(x=as.matrix(Real_Klein))

scran_012<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta = 1,inputBeta=scransf/mean(scransf)*0.12,inputtrim = 0.01)
scran_006<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta =1,inputBeta=scransf/mean(scransf)*0.06,inputtrim = 0.01)
scran_003<-trytt2(counts=as.matrix(Real_Klein),inputMeanBeta = 1,inputBeta=scransf/mean(scransf)*0.03,inputtrim = 0.01)




real=SingleCellExperiment(assays=list(counts=as.matrix(Real_Klein)))

sce_trim001_012=SingleCellExperiment(assays=list(counts=trim001_012$SCE@assays@.xData$data$counts))
sce_trim001_006=SingleCellExperiment(assays=list(counts=trim001_006$SCE@assays@.xData$data$counts))
sce_trim001_003=SingleCellExperiment(assays=list(counts=trim001_003$SCE@assays@.xData$data$counts))

sce_scran_012=SingleCellExperiment(assays=list(counts=scran_012$SCE@assays@.xData$data$counts))
sce_scran_006=SingleCellExperiment(assays=list(counts=scran_006$SCE@assays@.xData$data$counts))
sce_scran_003=SingleCellExperiment(assays=list(counts=scran_003$SCE@assays@.xData$data$counts))




trim001_012$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_006$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_003$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_012$SCE@assays@.xData$data$BCV<-NULL
trim001_006$SCE@assays@.xData$data$BCV<-NULL
trim001_003$SCE@assays@.xData$data$BCV<-NULL
trim001_012$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_003$SCE@assays@.xData$data$Single_Genemean_mat<-NULL


scran_012$SCE@assays@.xData$data$TrueCounts<-NULL
scran_006$SCE@assays@.xData$data$TrueCounts<-NULL
scran_003$SCE@assays@.xData$data$TrueCounts<-NULL
scran_012$SCE@assays@.xData$data$BCV<-NULL
scran_006$SCE@assays@.xData$data$BCV<-NULL
scran_003$SCE@assays@.xData$data$BCV<-NULL
scran_012$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_003$SCE@assays@.xData$data$Single_Genemean_mat<-NULL

trim001_006$Est_params@counts.origin<-matrix()
trim001_006$Est_params@counts.norm.TC<-matrix()
trim001_012$Est_params@counts.origin<-matrix()
trim001_012$Est_params@counts.norm.TC<-matrix()
trim001_003$Est_params@counts.origin<-matrix()
trim001_003$Est_params@counts.norm.TC<-matrix()

scran_006$Est_params@counts.origin<-matrix()
scran_006$Est_params@counts.norm.TC<-matrix()
scran_012$Est_params@counts.origin<-matrix()
scran_012$Est_params@counts.norm.TC<-matrix()
scran_003$Est_params@counts.origin<-matrix()
scran_003$Est_params@counts.norm.TC<-matrix()



SCElist_Klein<-list(
    real=real,
    trim001_003=sce_trim001_003,trim001_006=sce_trim001_006,trim001_012=sce_trim001_012,
    scran_003=sce_scran_003,scran_006=sce_scran_006,scran_012=sce_scran_012
)



library(SummarizedExperiment)
library(SingleCellExperiment)
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
names(SCElist_Klein)<-c('Real data',"trim1%_3%", "trim1%_6%", "trim1%_12%","scran_3%" ,  "scran_6%"  , "scran_12%")
library(ggplot2)
source("E:/RNAseqProject/PROJECT_scRNAseq/FIGURE_DROPOUT/DROPOUT_MAIN/DROPOUT_FUN.r")
meanvar<-plot_MEANVAR_v2(sces=SCElist_Klein,MAIN='',CAPTION='Klein case study',legend.position=c(0.5,0.9))
meanvar_diff<-plot_MEANVAR_diff_v2(sces=SCElist_Klein,MAIN='',CAPTION='Klein case study')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    meanvar$mean.zeros,
    meanvar_diff$mean.zeros,
    
    meanvar$mean.var,
    meanvar_diff$mean.var,
    
    meanvar$z.gene,
    meanvar_diff$z.gene,
    
    meanvar$z.cell,
    meanvar_diff$z.cell,
    
    ncol=2,nrow=4,common.legend = TRUE, legend="bottom")

qq

ggsave(file="E:/RNAseqProject/SIMULATION/SIM_EXPLORE/Klein_simexplore.pdf",device='pdf',plot=qq,width=8,height=10)



#save.image(file="E:/RNAseqProject/Klein2015/Klein_simexplore/Klein_simexplore.RData")
load("E:/RNAseqProject/SIMULATION/SIM_EXPLORE/Klein_simexplore.RData")
rm(meanvar)
rm(meanvar_diff)
rm(bay_Klein)
rm(qq)
rm(scran_003)
rm(scran_006)
rm(scran_012)
rm(trim001_003)
rm(trim001_006)
rm(trim001_012)


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    # meanvar$mean.zeros,
    # meanvar_diff$mean.zeros,
    # 
    # meanvar$mean.var,
    # meanvar_diff$mean.var,
    
    #meanvar$z.gene,
    #meanvar_diff$z.gene,
    
    as_ggplot(meanvar$z.cell),
    meanvar_diff$z.cell,
    
    ncol=2,nrow=1,common.legend = TRUE, legend="top")

qq
ggsave(file="E:/RNAseqProject/Illustrator_bayNorm/newFIG/sub_Klein.pdf",device='pdf',plot=qq,width=8.2,height=3.5)

########Tung N1########################
load('Tung_N1_SIM')


library(ggplot2)

source("E:/RNAseqProject/MANY_SAVE_PATH.R")
source("E:/RNAseqProject/PROJECT_scRNAseq/FIGURE_DROPOUT/DROPOUT_MAIN/DROPOUT_FUN.r")

library(SingleCellExperiment)

real=SingleCellExperiment(assays=list(counts=as.matrix(N1_DAT)))

spikein_012<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 1,inputBeta=bayNorm_mean_N1$BETA)
spikein_006<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 1,inputBeta=bayNorm_mean_N1$BETA/mean(bayNorm_mean_N1$BETA)*0.06)
spikein_024<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 1,inputBeta=bayNorm_mean_N1$BETA/mean(bayNorm_mean_N1$BETA)*0.24)

trim001_012<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 0.12,inputBeta=NULL,inputtrim = 0.01)
trim001_006<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 0.06,inputBeta=NULL,inputtrim = 0.01)
trim001_024<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 0.24,inputBeta=NULL,inputtrim = 0.01)

library(scran)
scransf<-computeSumFactors(x=N1_DAT)

scran_012<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 1,inputBeta=scransf/mean(scransf)*0.12,inputtrim = 0.01)
scran_006<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta =1,inputBeta=scransf/mean(scransf)*0.06,inputtrim = 0.01)
scran_024<-trytt2(counts=as.matrix(N1_DAT),inputMeanBeta = 1,inputBeta=scransf/mean(scransf)*0.24,inputtrim = 0.01)




real=SingleCellExperiment(assays=list(counts=as.matrix(N1_DAT)))
sce_spikein_012=SingleCellExperiment(assays=list(counts=spikein_012$SCE@assays@.xData$data$counts))
sce_spikein_006=SingleCellExperiment(assays=list(counts=spikein_006$SCE@assays@.xData$data$counts))
sce_spikein_024=SingleCellExperiment(assays=list(counts=spikein_024$SCE@assays@.xData$data$counts))

sce_trim001_012=SingleCellExperiment(assays=list(counts=trim001_012$SCE@assays@.xData$data$counts))
sce_trim001_006=SingleCellExperiment(assays=list(counts=trim001_006$SCE@assays@.xData$data$counts))
sce_trim001_024=SingleCellExperiment(assays=list(counts=trim001_024$SCE@assays@.xData$data$counts))

sce_scran_012=SingleCellExperiment(assays=list(counts=scran_012$SCE@assays@.xData$data$counts))
sce_scran_006=SingleCellExperiment(assays=list(counts=scran_006$SCE@assays@.xData$data$counts))
sce_scran_024=SingleCellExperiment(assays=list(counts=scran_024$SCE@assays@.xData$data$counts))






spikein_012$SCE@assays@.xData$data$TrueCounts<-NULL
spikein_006$SCE@assays@.xData$data$TrueCounts<-NULL
spikein_024$SCE@assays@.xData$data$TrueCounts<-NULL
spikein_012$SCE@assays@.xData$data$BCV<-NULL
spikein_006$SCE@assays@.xData$data$BCV<-NULL
spikein_024$SCE@assays@.xData$data$BCV<-NULL
spikein_012$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
spikein_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
spikein_024$SCE@assays@.xData$data$Single_Genemean_mat<-NULL

trim001_012$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_006$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_024$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_012$SCE@assays@.xData$data$BCV<-NULL
trim001_006$SCE@assays@.xData$data$BCV<-NULL
trim001_024$SCE@assays@.xData$data$BCV<-NULL
trim001_012$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_024$SCE@assays@.xData$data$Single_Genemean_mat<-NULL


scran_012$SCE@assays@.xData$data$TrueCounts<-NULL
scran_006$SCE@assays@.xData$data$TrueCounts<-NULL
scran_024$SCE@assays@.xData$data$TrueCounts<-NULL
scran_012$SCE@assays@.xData$data$BCV<-NULL
scran_006$SCE@assays@.xData$data$BCV<-NULL
scran_024$SCE@assays@.xData$data$BCV<-NULL
scran_012$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_024$SCE@assays@.xData$data$Single_Genemean_mat<-NULL

spikein_006$Est_params@counts.origin<-matrix()
spikein_006$Est_params@counts.norm.TC<-matrix()
spikein_012$Est_params@counts.origin<-matrix()
spikein_012$Est_params@counts.norm.TC<-matrix()
spikein_024$Est_params@counts.origin<-matrix()
spikein_024$Est_params@counts.norm.TC<-matrix()

trim001_006$Est_params@counts.origin<-matrix()
trim001_006$Est_params@counts.norm.TC<-matrix()
trim001_012$Est_params@counts.origin<-matrix()
trim001_012$Est_params@counts.norm.TC<-matrix()
trim001_024$Est_params@counts.origin<-matrix()
trim001_024$Est_params@counts.norm.TC<-matrix()

scran_006$Est_params@counts.origin<-matrix()
scran_006$Est_params@counts.norm.TC<-matrix()
scran_012$Est_params@counts.origin<-matrix()
scran_012$Est_params@counts.norm.TC<-matrix()
scran_024$Est_params@counts.origin<-matrix()
scran_024$Est_params@counts.norm.TC<-matrix()



SCElist_TungN1<-list(
    real=real,
    spikein_006=sce_spikein_006,spikein_012=sce_spikein_012,spikein_024=sce_spikein_024,
    trim001_006=sce_trim001_006,trim001_012=sce_trim001_012,trim001_024=sce_trim001_024,
    scran_006=sce_scran_006,scran_012=sce_scran_012,scran_024=sce_scran_024
)


library(SummarizedExperiment)
library(SingleCellExperiment)
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

names(SCElist_TungN1)<-c('Real data',"spikein_6%", "spikein_12%", "spikein_24%","trim1%_6%", "trim1%_12%", "trim1%_24%", "scran_6%" ,"scran_12%"  , "scran_24%")

source("DROPOUT_FUN.r")
meanvar<-plot_MEANVAR_v2(sces=SCElist_TungN1,MAIN='',CAPTION='Tung case study (NA19098)',legend.position=c(0.5,0.9))
meanvar_diff<-plot_MEANVAR_diff_v2(sces=SCElist_TungN1,MAIN='',CAPTION='Tung case study (NA19098)')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    # meanvar$mean.zeros,
    # meanvar_diff$mean.zeros,
    # 
    # meanvar$mean.var,
    # meanvar_diff$mean.var,
    # 
    # meanvar$z.gene,
    # meanvar_diff$z.gene,
    
    as_ggplot(meanvar$z.cell),
    meanvar_diff$z.cell,
    
    ncol=2,nrow=1,common.legend = TRUE, legend="bottom")

qq




########Torre########################
load('Torre_sim.RData')
library(ggplot2)

source("E:/RNAseqProject/MANY_SAVE_PATH.R")
source("E:/RNAseqProject/PROJECT_scRNAseq/FIGURE_DROPOUT/DROPOUT_MAIN/DROPOUT_FUN.r")
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
source("E:/RNAseqProject/MANY_SAVE_PATH.R")
source("E:/RNAseqProject/PROJECT_scRNAseq/FIGURE_DROPOUT/DROPOUT_MAIN/DROPOUT_FUN.r")



real=SingleCellExperiment(assays=list(counts=as.matrix(Torre_drop_sub)))


GAPDH_003<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 1,inputBeta=bay_out$BETA)
GAPDH_006<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 1,inputBeta=bay_out$BETA/mean(bay_out$BETA)*0.06)
GAPDH_0015<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 1,inputBeta=bay_out$BETA/mean(bay_out$BETA)*0.015)

trim001_003<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 0.03,inputBeta=NULL,inputtrim = 0.01)
trim001_006<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 0.06,inputBeta=NULL,inputtrim = 0.01)
trim001_0015<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 0.015,inputBeta=NULL,inputtrim = 0.01)

library(scran)
scransf<-computeSumFactors(x=Torre_drop_sub)

b1<-scransf/mean(scransf)*0.03
summary(b1)
b1[b1>=1]=max(b1[b1<1])
scran_003<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 1,inputBeta=b1,inputtrim = 0.01)

b2<-scransf/mean(scransf)*0.06
summary(b2)
b2[b2>=1]=max(b2[b2<1])
scran_006<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta =1,inputBeta=b2,inputtrim = 0.01)

b3<-scransf/mean(scransf)*0.015
summary(b3)
b3[b3>=1]=max(b3[b3<1])
scran_0015<-trytt2(counts=as.matrix(Torre_drop_sub),inputMeanBeta = 1,inputBeta=b3,inputtrim = 0.01)


real=SingleCellExperiment(assays=list(counts=as.matrix(Torre_drop_sub)))


sce_GAPDH_003=SingleCellExperiment(assays=list(counts=GAPDH_003$SCE@assays@.xData$data$counts))
sce_GAPDH_006=SingleCellExperiment(assays=list(counts=GAPDH_006$SCE@assays@.xData$data$counts))
sce_GAPDH_0015=SingleCellExperiment(assays=list(counts=GAPDH_0015$SCE@assays@.xData$data$counts))


sce_trim001_003=SingleCellExperiment(assays=list(counts=trim001_003$SCE@assays@.xData$data$counts))
sce_trim001_006=SingleCellExperiment(assays=list(counts=trim001_006$SCE@assays@.xData$data$counts))
sce_trim001_0015=SingleCellExperiment(assays=list(counts=trim001_0015$SCE@assays@.xData$data$counts))

sce_scran_003=SingleCellExperiment(assays=list(counts=scran_003$SCE@assays@.xData$data$counts))
sce_scran_006=SingleCellExperiment(assays=list(counts=scran_006$SCE@assays@.xData$data$counts))
sce_scran_0015=SingleCellExperiment(assays=list(counts=scran_0015$SCE@assays@.xData$data$counts))

GAPDH_003$SCE@assays@.xData$data$TrueCounts<-NULL
GAPDH_006$SCE@assays@.xData$data$TrueCounts<-NULL
GAPDH_0015$SCE@assays@.xData$data$TrueCounts<-NULL
GAPDH_003$SCE@assays@.xData$data$BCV<-NULL
GAPDH_006$SCE@assays@.xData$data$BCV<-NULL
GAPDH_0015$SCE@assays@.xData$data$BCV<-NULL
GAPDH_003$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
GAPDH_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
GAPDH_0015$SCE@assays@.xData$data$Single_Genemean_mat<-NULL


trim001_003$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_006$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_0015$SCE@assays@.xData$data$TrueCounts<-NULL
trim001_003$SCE@assays@.xData$data$BCV<-NULL
trim001_006$SCE@assays@.xData$data$BCV<-NULL
trim001_0015$SCE@assays@.xData$data$BCV<-NULL
trim001_003$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
trim001_0015$SCE@assays@.xData$data$Single_Genemean_mat<-NULL


scran_003$SCE@assays@.xData$data$TrueCounts<-NULL
scran_006$SCE@assays@.xData$data$TrueCounts<-NULL
scran_0015$SCE@assays@.xData$data$TrueCounts<-NULL
scran_003$SCE@assays@.xData$data$BCV<-NULL
scran_006$SCE@assays@.xData$data$BCV<-NULL
scran_0015$SCE@assays@.xData$data$BCV<-NULL
scran_003$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_006$SCE@assays@.xData$data$Single_Genemean_mat<-NULL
scran_0015$SCE@assays@.xData$data$Single_Genemean_mat<-NULL

GAPDH_006$Est_params@counts.origin<-matrix()
GAPDH_006$Est_params@counts.norm.TC<-matrix()
GAPDH_003$Est_params@counts.origin<-matrix()
GAPDH_003$Est_params@counts.norm.TC<-matrix()
GAPDH_0015$Est_params@counts.origin<-matrix()
GAPDH_0015$Est_params@counts.norm.TC<-matrix()

trim001_006$Est_params@counts.origin<-matrix()
trim001_006$Est_params@counts.norm.TC<-matrix()
trim001_003$Est_params@counts.origin<-matrix()
trim001_003$Est_params@counts.norm.TC<-matrix()
trim001_0015$Est_params@counts.origin<-matrix()
trim001_0015$Est_params@counts.norm.TC<-matrix()

scran_006$Est_params@counts.origin<-matrix()
scran_006$Est_params@counts.norm.TC<-matrix()
scran_003$Est_params@counts.origin<-matrix()
scran_003$Est_params@counts.norm.TC<-matrix()
scran_0015$Est_params@counts.origin<-matrix()
scran_0015$Est_params@counts.norm.TC<-matrix()



SCElist_Torre<-list(
    real=real,
    GAPDH_0015=sce_GAPDH_0015,GAPDH_003=sce_GAPDH_003,GAPDH_006=sce_GAPDH_006,
    trim001_0015=sce_trim001_0015,trim001_003=sce_trim001_003,trim001_006=sce_trim001_006,
    scran_0015=sce_scran_0015,scran_003=sce_scran_003,scran_006=sce_scran_006
)


names(SCElist_Torre)<-c("Real data" ,  "GAPDH_1.5%",   "GAPDH_3%",    "GAPDH_6%","trim1%_1.5%", "trim1%_3%" , "trim1%_6%",  "scran_1.5%","scran_3%","scran_6%" )



library(SummarizedExperiment)
library(SingleCellExperiment)
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
meanvar<-plot_MEANVAR_v2(sces=SCElist_Torre,MAIN='',CAPTION='Torre case study',legend.position='top')
meanvar_diff<-plot_MEANVAR_diff_v2(sces=SCElist_Torre,MAIN='',CAPTION='Torre case study')


library(gridExtra)
library(ggpubr)
qq<-ggarrange(
    # meanvar$mean.zeros,
    # meanvar_diff$mean.zeros,
    # 
    # meanvar$mean.var,
    # meanvar_diff$mean.var,
    # 
    # meanvar$z.gene,
    # meanvar_diff$z.gene,
    
    as_ggplot(meanvar$z.cell),
    meanvar_diff$z.cell,
    
    ncol=2,nrow=1,common.legend = TRUE, legend="bottom")

qq
