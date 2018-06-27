load("Grun_2i_norms.RData")
load("Grun_serum_norms.RData")
load("smFISH_norm_load.RData")



SCnorm_2i<-SCnorm_2i_ori@metadata$NormalizedData
SCnorm_serum<-SCnorm_serum_ori@metadata$NormalizedData


MAGIC_2i<-MAGIC_2i/mean(bayNorm_SC_2i$BETA)
DCA_2i<-DCA_2i/mean(bayNorm_SC_2i$BETA)
scImpute_SC_2i<-scImpute_SC_2i/mean(bayNorm_SC_2i$BETA)
SAVER_SC_2i_array<-SAVER_SC_2i_array/mean(bayNorm_SC_2i$BETA)
SCnorm_2i<-SCnorm_2i/mean(bayNorm_SC_2i$BETA)
MAGIC_2i<-MAGIC_2i/mean(bayNorm_SC_2i$BETA)

DCA_serum<-DCA_serum/mean(bayNorm_SC_serum$BETA)
scImpute_SC_serum<-scImpute_SC_serum/mean(bayNorm_SC_serum$BETA)
SAVER_SC_serum_array<-SAVER_SC_serum_array/mean(bayNorm_SC_serum$BETA)
SCnorm_serum<-SCnorm_serum/mean(bayNorm_SC_serum$BETA)
MAGIC_serum<-MAGIC_serum/mean(bayNorm_SC_serum$BETA)


CV_MAGIC_2i<-apply(MAGIC_2i,1,sd)/rowMeans(MAGIC_2i)
CV_MAGIC_serum<-apply(MAGIC_serum,1,sd)/rowMeans(MAGIC_serum)
CV_DCA_2i<-apply(DCA_2i,1,sd)/rowMeans(DCA_2i)
CV_DCA_serum<-apply(DCA_serum,1,sd)/rowMeans(DCA_serum)


set.seed(1230000)
CV_fun <- function(x) sqrt(var(x))/mean(x)

#smFISH bootstrapping: 2i#######
library(foreach)
boot_2imat_smFISH<-matrix(nrow=1000,ncol=9)
colnames(boot_2imat_smFISH)<-Fig3Gene
for(i in 1:length(Fig3Gene)){
  print(i)
  if(Fig3Gene[i]!='Pou5f1'){
    
    for(j in 1:length(Grun_2i_norm_list)){
      
      if(length(intersect(Fig3Gene[i],colnames(Grun_2i_norm_list[[j]])))==1){
        asd<-which(colnames(Grun_2i_norm_list[[j]])==Fig3Gene[i])
        boot_2imat_smFISH[,i]<-foreach(qq=1:1000,.combine=c)%do%{
          qq_temp<-CV_fun(sample(Grun_2i_norm_list[[j]][,asd],replace=T))
          return(qq_temp)
        }
      }
    }
    
  }else if(Fig3Gene[i]=='Pou5f1'){
    #useeeee<-c(Grun_2i_norm_list$Notch1_Pou5f1_2i[,2],Grun_2i_norm_list$Pou5f1_Klf4_2i[,1],Grun_2i_norm_list$Sox2_Pou5f1_2i[,2])
    useeeee<-c(Grun_2i_norm_list$Pou5f1_Klf4_2i[,1])
    
    boot_2imat_smFISH[,9]<-foreach(qq=1:1000,.combine=c)%do%{
      qq_temp<-CV_fun(sample(useeeee,replace=T))
      return(qq_temp)
    }
  }
}

#smFISH bootstrapping: serum######
library(foreach)
boot_serummat_smFISH<-matrix(nrow=1000,ncol=9)
colnames(boot_serummat_smFISH)<-Fig3Gene

for(i in 1:length(Fig3Gene)){
  print(i)
  if(Fig3Gene[i]!='Pou5f1'){
    
    for(j in 1:length(Grun_S_norm_list)){
      
      if(length(intersect(Fig3Gene[i],colnames(Grun_S_norm_list[[j]])))==1){
        asd<-which(colnames(Grun_S_norm_list[[j]])==Fig3Gene[i])
        boot_serummat_smFISH[,i]<-foreach(qq=1:1000,.combine=c)%do%{
          qq_temp<-CV_fun(sample(Grun_S_norm_list[[j]][,asd],replace=T))
          return(qq_temp)
        }
      }
    }
    
  }else if(Fig3Gene[i]=='Pou5f1'){
    #useeeee<-c(Grun_S_norm_list$Pou5f1_Klf4_J1S[,1],Grun_S_norm_list$Pou5f1_Notch1_S[,1])
    useeeee<-c(Grun_S_norm_list$Pou5f1_Klf4_J1S[,1])
    boot_serummat_smFISH[,9]<-foreach(qq=1:1000,.combine=c)%do%{
      qq_temp<-CV_fun(sample(useeeee,replace=T))
      return(qq_temp)
    }
  }
}
sum(is.na(boot_serummat_smFISH))

###bayNorm 2i:########
boot_2imat_bayNorm<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(bayNorm_SC_2i$Bay_array)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(bayNorm_SC_2i$Bay_array[qqind,,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_bayNorm)<-Fig3Gene

###bayNorm serum:####
boot_serummat_bayNorm<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(bayNorm_SC_serum$Bay_array)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(bayNorm_SC_serum$Bay_array[qqind,,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_bayNorm)<-Fig3Gene


######SCnorm ####

boot_2imat_scnorm<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(SCnorm_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(SCnorm_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_scnorm)<-Fig3Gene

boot_serummat_scnorm<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(SCnorm_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(SCnorm_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_scnorm)<-Fig3Gene


###scran########

boot_2imat_scran<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(scran_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(scran_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_scran)<-Fig3Gene

boot_serummat_scran<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(scran_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(scran_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_scran)<-Fig3Gene


#####RPM#########

boot_2imat_RPM<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(RPM_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(RPM_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_RPM)<-Fig3Gene

boot_serummat_RPM<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(RPM_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(RPM_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_RPM)<-Fig3Gene

######TMM#########
boot_2imat_TMM<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(TMM_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(TMM_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_TMM)<-Fig3Gene

boot_serummat_TMM<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(TMM_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(TMM_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_TMM)<-Fig3Gene


#######DESeq#######
boot_2imat_DESeq<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(DESeq_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(DESeq_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_DESeq)<-Fig3Gene

boot_serummat_DESeq<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(DESeq_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(DESeq_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_DESeq)<-Fig3Gene



#######SAVER#######
boot_2imat_SAVER<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(SAVER_SC_2i$estimate)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(as.vector(SAVER_SC_2i_array[qqind,,]),replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_SAVER)<-Fig3Gene

boot_serummat_SAVER<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(SAVER_SC_serum$estimate)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(as.vector(SAVER_SC_serum_array[qqind,,]),replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_SAVER)<-Fig3Gene



#######RAW_BETA#######
boot_2imat_RAW_BETA<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(RAW_BETA_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(RAW_BETA_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_RAW_BETA)<-Fig3Gene

boot_serummat_RAW_BETA<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(RAW_BETA_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(RAW_BETA_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_RAW_BETA)<-Fig3Gene

#######scImpute#######
boot_2imat_scImpute<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(scImpute_SC_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(scImpute_SC_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_scImpute)<-Fig3Gene

boot_serummat_scImpute<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(scImpute_SC_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(scImpute_SC_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_scImpute)<-Fig3Gene


#######MAGIC#######
boot_2imat_MAGIC<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(MAGIC_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(MAGIC_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_MAGIC)<-Fig3Gene

boot_serummat_MAGIC<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(MAGIC_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(MAGIC_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_MAGIC)<-Fig3Gene

#######DCA#######
boot_2imat_DCA<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(DCA_2i)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(DCA_2i[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_2imat_DCA)<-Fig3Gene

boot_serummat_DCA<-foreach(i=1:9,.combine=cbind)%:%
  foreach(j=1:1000,.combine=c)%do%{
    qqind<-which(rownames(DCA_serum)==Fig3Gene[i])
    qqtemp<-CV_fun(sample(DCA_serum[qqind,],replace=T))
    return(qqtemp)
  }
colnames(boot_serummat_DCA)<-Fig3Gene

#Store into list#####
STORE_BOOTCV_list<-list(24)
qqq<-ls()
qqq2<-grep(ls(),pattern='boot_')
length(qqq2)

for(i in 1:24){
  STORE_BOOTCV_list[[i]]<-get(qqq[qqq2[i]])
}

names(STORE_BOOTCV_list)<-qqq[qqq2]


#new MAGIC, default SAVER, DCA
save(STORE_BOOTCV_list,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_smFISH_meanBETA/Bootstrapping_CV_V3tr_default_divmeanbeta.RData")
