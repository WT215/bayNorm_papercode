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


save.image('Klein_bayNorm.RData')