#Prepare smFISH data######
Torre_FISH<-read.table(file="fishSubset.txt")


library(foreach)
FISH_mu<-foreach(i=1:26,.combine=c)%do%{
    x<-Torre_FISH[,i]
    qq<-mean(x[!is.na(x)])
}
names(FISH_mu)<-colnames(Torre_FISH)

FISH_cv<-foreach(i=1:26,.combine=c)%do%{
    x<-Torre_FISH[,i]
    qq<-mean(x[!is.na(x)])
    sd<-sd(x[!is.na(x)])
    cv<-sd/qq
    return(cv)
}
names(FISH_cv)<-colnames(Torre_FISH)

summary(Torre_FISH[,'GAPDH'])
GAPDH_10QUANTILE<-quantile(Torre_FISH[,'GAPDH'],probs=c(0.1,0.9))
celldrop<-c(which(Torre_FISH[,'GAPDH']<GAPDH_10QUANTILE[1]) , which(Torre_FISH[,'GAPDH']>GAPDH_10QUANTILE[2]))

Torre_FISH_sub<-Torre_FISH[-celldrop,]
GAPDH_sf<-Torre_FISH_sub[,'GAPDH']/mean(Torre_FISH_sub[,'GAPDH'])

Torre_FISH_sub_norm<-Torre_FISH_sub/GAPDH_sf

sum(is.na(Torre_FISH_sub_norm))


FISH_MEAN_norm<-apply(Torre_FISH_sub_norm,2,function(x){mean(x[!is.na(x)])})

#Saver preprocessed smFISH data for further studying
save.image("smFISH.RData")


#Preprocess drop seq data####
#Functions and data were kindly provided by Mo Huang
df.upm <- as.matrix(read.table("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/Data/GSE99330_dropseqUPM.txt", row.names = 1, header = TRUE))
Torre_dropseq <- sweep(df.upm, 2, apply(df.upm, 2, function(x) min(x[x!= 0])), "/")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

corverage<-apply(Torre_dropseq,2,function(x){length(which(x!=0))})
GAPDH_raw<-Torre_dropseq['GAPDH',]

Torre_drop_sub<-Torre_dropseq[,intersect(which(corverage>2000),which(GAPDH_raw!=0))]
Torre_drop_sub<-Torre_drop_sub[-which(rowMeans(Torre_drop_sub)<0.1),]


(smFISH_genename<-intersect(rownames(Torre_drop_sub),names(smFISH_list)))



(xx<-(unlist(lapply(smFISH_list,mean))[smFISH_genename]))
(yy<-rowMeans(Torre_drop_sub[smFISH_genename,]))
plot(xx,yy,log='xy')
text(xx,yy,names(xx),pos=1)
dropgene<-c('VGF','MITF','SOX10')
dropind<-which(smFISH_genename %in% dropgene)
smFISH_genename
(xx<-(unlist(lapply(smFISH_list,mean))[smFISH_genename])[-dropind])
(yy<-rowMeans(Torre_drop_sub[smFISH_genename,])[-dropind])
lmm<-lm(yy~xx)
lmm_summa<-summary(lmm)


#Sup Fig 
plot(xx,yy,xlab='Mean expression of FISH counts',ylab='Mean expression of raw data',main='Data from Torre et al study',pch=16,cex.main=0.95)
abline(lmm,lty=2)
#text(xx,yy,names(xx),pos=1)
legend('topleft',legend=c(paste('Estimated beta=',round(coef(lmm)[2],4)),paste('adj R2=',round(lmm_summa$adj.r.squared,4)),paste('cor=',round(cor(xx,yy),4))),bty='n')

(lmm_summa<-summary(lmm))
(MEANBETA<-coef(lmm)[2])


BETA<-Torre_drop_sub['GAPDH',]/median(Torre_drop_sub['GAPDH',])*MEANBETA
summary(BETA)

BETA[BETA<=0]<-min(BETA[BETA>0])
BETA[BETA>=1]<-max(BETA[BETA<1])
summary(BETA)

Torre_drop_sub<-Torre_drop_sub[-which(rownames(Torre_drop_sub)=='GAPDH'),]

#Then we feed Torre_drop_sub into different normalization methods
save.image('Load_Torre.RData')



