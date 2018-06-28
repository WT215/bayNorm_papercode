library(fdrtool)
library(ggplot2)
#########################################
########Vahid's T_stat##################
#########################################

#Modified by Wenhao on 14/08/2017. Use p.adjust instead of fdrtool:
#Return additional results, both came from fdrtool:
#1. TFDR
#2. PValue


#Modified by Wenhao on 09/08/2017. Use p.adjust instead of fdrtool
#Return:
#1. TSTAT
#2. PVAL_raw: raw p-value of each test
#3. PVAL_fdr: adjusted p-value of each test


#home=Sys.getenv("HOME")

#source(paste(home,"/ANALYSIS/BAYESIAN_NORM/functions_1408.R",sep=''))

tStatAnalysis=function(cells,ctrls, list_mode=NULL, verbose=F, ...)
{
 #source("~/Box Sync/UMIs/Functions/functions_1408.R")


  if(length(dim(cells))==2 & length(dim(ctrls))==2){
    cells<-array(cells,dim=c(dim(cells),1))
    ctrls<-array(ctrls,dim=c(dim(ctrls),1))
  }


##format input##
 if(is.null(list_mode)==T)
 {
  CELLS=cells
  CTRLS=ctrls
 }
 if(is.null(list_mode)==F)
 {
  if((is.null(dim(list_mode))==T)&&(length(list_mode) > 0))
  {
   CELLS=cells[list_mode,,]
   CTRLS=ctrls[list_mode,,]
  }
  else if(class(list_mode)=="list")
  {
   CELLS=array(dim=c(length(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(names(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
   CTRLS=array(dim=c(length(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(names(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))
   for(i in 1:lenght(list_mode))
   {
    for(j in 1:dim(cells)[3])
    {
     GOCAT=list_mode[[i]]
     GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
     if(length(GOCAT) > 1)
     {
      CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
      CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
     }
     else if(length(GOCAT)==1)
     {
      CELLS[i,,j]=cells[GOCAT,,j]
      CTRLS[i,,j]=ctrls[GOCAT,,j]
     }
     else
     {
      CELLS[i,,j]=rep(0,dim(cells)[2])
      CTRLS[i,,j]=rep(0,dim(ctrls)[2])
     }
    }
   }
  }
  else if((length(dim(list_mode)==2))&&((list_mode==1)||(list_mode==0)))
  {
   CELLS=array(dim=c(ncol(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
   CTRLS=array(dim=c(ncol(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))

   for(i in 1:ncol(list_mode))
   {
    for(j in 1:dim(cells)[3])
    {
     GOCAT=row.names(list_mode)[which(list_mode[,i]==1)]
     GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
     if(length(GOCAT) > 1)
     {
      CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
      CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
     }
     else if(length(GOCAT)==1)
     {
      CELLS[i,,j]=cells[GOCAT,,j]
      CTRLS[i,,j]=ctrls[GOCAT,,j]
     }
     else
     {
      CELLS[i,,j]=rep(0,dim(cells)[2])
      CTRLS[i,,j]=rep(0,dim(ctrls)[2])
     }
    }
   }
  }
  else {stop("Unrecognise list_mode format")}
 }
##Create output tables##
 TSTAT=matrix(NA,dim(CELLS)[1],dim(CELLS)[2])
 PVAL_raw=matrix(NA,dim(CELLS)[1],dim(CELLS)[2])
 PVAL_fdr=matrix(NA,dim(CELLS)[1],dim(CELLS)[2])

 TFDR=matrix(NA,dim(CELLS)[1], dim(CELLS)[2])
 Tfdr=matrix(NA,dim(CELLS)[1], dim(CELLS)[2])
 PValue=matrix(NA,dim(CELLS)[1], dim(CELLS)[2])
 Param=matrix(NA,6, dim(CELLS)[2])

 row.names(TSTAT)=dimnames(CELLS)[[1]]
 colnames(TSTAT)=dimnames(CELLS)[[2]]

 row.names(PVAL_raw)=dimnames(CELLS)[[1]]
 colnames(PVAL_raw)=dimnames(CELLS)[[2]]

 row.names(PVAL_fdr)=dimnames(CELLS)[[1]]
 colnames(PVAL_fdr)=dimnames(CELLS)[[2]]

 row.names(TFDR)=dimnames(CELLS)[[1]]
 colnames(TFDR)=dimnames(CELLS)[[2]]
 row.names(Tfdr)=dimnames(CELLS)[[1]]
 colnames(Tfdr)=dimnames(CELLS)[[2]]
 row.names(PValue)=dimnames(CELLS)[[1]]
 colnames(PValue)=dimnames(CELLS)[[2]]
 colnames(Param)=dimnames(CELLS)[[2]]

##Compute T-statistics ##
 m = dim(CELLS)[3]
 n = dim(CTRLS)[2]*dim(CTRLS)[3]


 for(i in 1:ncol(CELLS))
 {
  if(verbose){print(i)}
  for (j in 1:dim(CELLS)[1]) {
  wtest =wilcox.test(CELLS[j,i,],CTRLS[j,,])
  TSTAT[j,i] = (wtest$statistic - m*n/2)/m/n*2
  PVAL_raw[j,i]<-wtest$p.value
  }
 }
 if(verbose){cat("\nDone with T statistics\n")}

 ##Compute FDR using fdrtool package ##
 library("fdrtool")

 for(i in 1:ncol(CELLS)){
   if(verbose){print(i)}
   fdr =  fdrtool(TSTAT[,i], statistic = "normal", plot = FALSE, verbose=FALSE)
   PVAL_fdr[,i]<-p.adjust(PVAL_raw[,i],method='fdr')
   TFDR[,i] = fdr$qval
   #Tfdr[,i] = fdr$lfdr
   PValue[,i] = fdr$pval
   #Param[,i] = fdr$param
 }
 #row.names(Param)=colnames(fdr$param)

 if(verbose){cat("\nDone with FDR calculation\n")}

##Format and return output talbes##
 #return(list(TSTAT,TFDR, Tfdr, PValue, Param))
 return(list(TSTAT=TSTAT,TFDR=TFDR,PValue=PValue,PVAL_raw=PVAL_raw,PVAL_fdr=PVAL_fdr))
}

#####
#####

tStatAnalysis_2groups =function(cells,ctrls, list_mode=NULL, verbose=F, plotFDR = FALSE, ...)
{


  if(length(dim(cells))==2 & length(dim(ctrls))==2){
    cells<-array(cells,dim=c(dim(cells),1))
    ctrls<-array(ctrls,dim=c(dim(ctrls),1))
  }

  #source("~/Box Sync/UMIs/Functions/functions_1408.R")

  ##format input##
  if(is.null(list_mode)==T)
  {
    CELLS=cells
    CTRLS=ctrls
  }
  if(is.null(list_mode)==F)
  {
    if((is.null(dim(list_mode))==T)&&(length(list_mode) > 0))
    {
      CELLS=cells[list_mode,,]
      CTRLS=ctrls[list_mode,,]
    }
    else if(class(list_mode)=="list")
    {
      CELLS=array(dim=c(length(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(names(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
      CTRLS=array(dim=c(length(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(names(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))
      for(i in 1:lenght(list_mode))
      {
        for(j in 1:dim(cells)[3])
        {
          GOCAT=list_mode[[i]]
          GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
          if(length(GOCAT) > 1)
          {
            CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
            CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
          }
          else if(length(GOCAT)==1)
          {
            CELLS[i,,j]=cells[GOCAT,,j]
            CTRLS[i,,j]=ctrls[GOCAT,,j]
          }
          else
          {
            CELLS[i,,j]=rep(0,dim(cells)[2])
            CTRLS[i,,j]=rep(0,dim(ctrls)[2])
          }
        }
      }
    }
    else if((length(dim(list_mode)==2))&&((list_mode==1)||(list_mode==0)))
    {
      CELLS=array(dim=c(ncol(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
      CTRLS=array(dim=c(ncol(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))

      for(i in 1:ncol(list_mode))
      {
        for(j in 1:dim(cells)[3])
        {
          GOCAT=row.names(list_mode)[which(list_mode[,i]==1)]
          GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
          if(length(GOCAT) > 1)
          {
            CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
            CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
          }
          else if(length(GOCAT)==1)
          {
            CELLS[i,,j]=cells[GOCAT,,j]
            CTRLS[i,,j]=ctrls[GOCAT,,j]
          }
          else
          {
            CELLS[i,,j]=rep(0,dim(cells)[2])
            CTRLS[i,,j]=rep(0,dim(ctrls)[2])
          }
        }
      }
    }
    else {stop("Unrecognise list_mode format")}
  }
  ##Create output tables##
  TSTAT=rep(NA,dim(CELLS)[1])
  PVAL_raw<-rep(NA,dim(CELLS)[1])
  PVAL_fdr<-rep(NA,dim(CELLS)[1])

  TFDR=rep(NA,dim(CELLS)[1])
  Tfdr=rep(NA,dim(CELLS)[1])
  PValue=rep(NA,dim(CELLS)[1])
  Param=rep(NA,6)

  names(TSTAT)=dimnames(CELLS)[[1]]

  names(PVAL_raw)<-dimnames(CELLS)[[1]]

  names(PVAL_fdr)<-dimnames(CELLS)[[1]]

  names(TFDR)=dimnames(CELLS)[[1]]

  names(Tfdr)=dimnames(CELLS)[[1]]

  names(PValue)=dimnames(CELLS)[[1]]



  ##Compute T-statistics ##
  m = dim(CELLS)[2]*dim(CELLS)[3]
  n = dim(CTRLS)[2]*dim(CTRLS)[3]

    for (j in 1:dim(CELLS)[1]) {
      wtest =wilcox.test(CELLS[j,,],CTRLS[j,,])
      TSTAT[j] = (wtest$statistic - m*n/2)/m/n*2
      PVAL_raw[j]<-wtest$p.value
    }

  if(verbose){cat("\nDone with T statistics\n")}

  ##Compute FDR using fdrtool package ##
  library("fdrtool")


    fdr =  fdrtool(TSTAT, statistic = "normal", plot = plotFDR, verbose=FALSE)
    TFDR = fdr$qval
    #Tfdr = fdr$lfdr
    PValue = fdr$pval
    #Param = fdr$param
  PVAL_fdr<-p.adjust(PVAL_raw,method='fdr')

#  row.names(Param)=colnames(fdr$param)

  if(verbose){cat("\nDone with FDR calculation\n")}

  ##Format and return output talbes##
  return(list(TSTAT=TSTAT,TFDR=TFDR,PValue=PValue,PVAL_raw=PVAL_raw,PVAL_fdr=PVAL_fdr))
}

#####
#####

FdrCurve<-function(Pval,TRUE_lab){
  Pval_s<-sort(Pval,decreasing=F)

  Fdddr<-foreach(i=1:length(Pval_s),.combine=rbind)%do%{
    S1<-which(Pval<Pval_s[i])
    S2<-which(TRUE_lab==1)
    fdrr<-(length(S1)-length(intersect(S1,S2)))/length(S1)
    hits<-length(S1)
    qq<-c(hits,fdrr)
    return(qq)
  }
  return(Fdddr)
}

####tStatAnalysis_2groups_2##########

tStatAnalysis_2groups_2 =function(cells,ctrls, list_mode=NULL, verbose=F, plotFDR = FALSE, ...)
{


  if(length(dim(cells))==2 & length(dim(ctrls))==2){
    cells<-array(cells,dim=c(dim(cells),1))
    ctrls<-array(ctrls,dim=c(dim(ctrls),1))
  }

  #source("~/Box Sync/UMIs/Functions/functions_1408.R")

  ##format input##
  if(is.null(list_mode)==T)
  {
    CELLS=cells
    CTRLS=ctrls
  }
  if(is.null(list_mode)==F)
  {
    if((is.null(dim(list_mode))==T)&&(length(list_mode) > 0))
    {
      CELLS=cells[list_mode,,]
      CTRLS=ctrls[list_mode,,]
    }
    else if(class(list_mode)=="list")
    {
      CELLS=array(dim=c(length(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(names(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
      CTRLS=array(dim=c(length(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(names(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))
      for(i in 1:lenght(list_mode))
      {
        for(j in 1:dim(cells)[3])
        {
          GOCAT=list_mode[[i]]
          GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
          if(length(GOCAT) > 1)
          {
            CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
            CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
          }
          else if(length(GOCAT)==1)
          {
            CELLS[i,,j]=cells[GOCAT,,j]
            CTRLS[i,,j]=ctrls[GOCAT,,j]
          }
          else
          {
            CELLS[i,,j]=rep(0,dim(cells)[2])
            CTRLS[i,,j]=rep(0,dim(ctrls)[2])
          }
        }
      }
    }
    else if((length(dim(list_mode)==2))&&((list_mode==1)||(list_mode==0)))
    {
      CELLS=array(dim=c(ncol(list_mode),ncol(cells),dim(cells)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(cells)[2],sep=''),paste("S",1:dim(cells)[3],sep='')))
      CTRLS=array(dim=c(ncol(list_mode),ncol(ctrls),dim(ctrls)[3]),dimnames=list(colnames(list_mode),paste("C",1:dim(ctrls)[2],sep=''),paste("S",1:dim(ctrls)[3],sep='')))

      for(i in 1:ncol(list_mode))
      {
        for(j in 1:dim(cells)[3])
        {
          GOCAT=row.names(list_mode)[which(list_mode[,i]==1)]
          GOCAT=GOCAT[which(GOCAT %in% dimnames(cells)[[1]])]
          if(length(GOCAT) > 1)
          {
            CELLS[i,,j]=colSums(cells[GOCAT,,j],na.rm=T)
            CTRLS[i,,j]=colSums(ctrls[GOCAT,,j],na.rm=T)
          }
          else if(length(GOCAT)==1)
          {
            CELLS[i,,j]=cells[GOCAT,,j]
            CTRLS[i,,j]=ctrls[GOCAT,,j]
          }
          else
          {
            CELLS[i,,j]=rep(0,dim(cells)[2])
            CTRLS[i,,j]=rep(0,dim(ctrls)[2])
          }
        }
      }
    }
    else {stop("Unrecognise list_mode format")}
  }
  ##Create output tables##
  TSTAT=rep(NA,dim(CELLS)[1])
  PVAL_raw<-rep(NA,dim(CELLS)[1])
  PVAL_fdr<-rep(NA,dim(CELLS)[1])

  TFDR=rep(NA,dim(CELLS)[1])
  Tfdr=rep(NA,dim(CELLS)[1])
  PValue=rep(NA,dim(CELLS)[1])
  Param=rep(NA,6)

  names(TSTAT)=dimnames(CELLS)[[1]]

  names(PVAL_raw)<-dimnames(CELLS)[[1]]

  names(PVAL_fdr)<-dimnames(CELLS)[[1]]

  names(TFDR)=dimnames(CELLS)[[1]]

  names(Tfdr)=dimnames(CELLS)[[1]]

  names(PValue)=dimnames(CELLS)[[1]]



  ##Compute T-statistics ##

  Result_list<-list()

  for(k in 1:dim(CELLS)[3]){
    print(paste('Sample',k,';',dim(CELLS)[3],'in total'))

    m = dim(CELLS)[2]
    n = dim(CTRLS)[2]


    for (j in 1:dim(CELLS)[1]) {
      wtest =wilcox.test(CELLS[j,,k],CTRLS[j,,k])
      TSTAT[j] = (wtest$statistic - m*n/2)/m/n*2
      PVAL_raw[j]<-wtest$p.value


    }
    fdr =  fdrtool(TSTAT, statistic = "normal", plot = plotFDR, verbose=FALSE)
    TFDR = fdr$qval
    PValue = fdr$pval
    PVAL_fdr<-p.adjust(PVAL_raw,method='fdr')

    Result_list[[k]]<-list(TSTAT=TSTAT,TFDR=TFDR,PValue=PValue,PVAL_raw=PVAL_raw,PVAL_fdr=PVAL_fdr)
  }



  if(verbose){cat("\nDone with T statistics\n")}

  ##Compute FDR using fdrtool package ##
  library("fdrtool")





  if(verbose){cat("\nDone with FDR calculation\n")}

  ##Format and return output talbes##
  return(Result_list=Result_list)
}



###wilcox_fun####
wilcox_fun<-function(cells,ctrls, list_mode=NULL, verbose=F, plotFDR = FALSE, ...){
    
   qq<- tStatAnalysis_2groups_2(cells,ctrls)
   
   qq2<-do.call(cbind,lapply(qq,function(x){return(x[[5]])}))
   qq2<-apply(qq2,1,median)
   names(qq2)<-rownames(cells)
   return(qq2)
}









    
wilcox_fun2<-function(cells,ctrls, list_mode=NULL, verbose=F, plotFDR = FALSE, ...){
    
    qq<- tStatAnalysis_2groups_2(cells,ctrls)
    
    qq1<-do.call(cbind,lapply(qq,function(x){return(x[[5]])}))
    qq2<-apply(qq1,1,median)
    names(qq2)<-rownames(cells)
    return(list(medianpval=qq2,matrixpval=qq1))
}    


wilcox_fun_array<-function(cells,ctrls){
    
    PVAL_raw<-rep(NA,dim(cells)[1])
    names(PVAL_raw)<-rownames(cells)
    

    for (j in 1:dim(cells)[1]) {
        wtest =wilcox.test(cells[j,,],ctrls[j,,])
        PVAL_raw[j]<-wtest$p.value
    }
    
    PVAL_fdr<-p.adjust(PVAL_raw,method='fdr')
    return(list( PVAL_fdr))
}  



kruskal_fun_array<-function(cells,ctrls){
    
    PVAL_raw<-rep(NA,dim(cells)[1])
    names(PVAL_raw)<-rownames(cells)
    
    for (j in 1:dim(cells)[1]) {
        ktest =kruskal.test(c(cells[j,,],ctrls[j,,]),g=c(rep(1,dim(cells)[2]*dim(cells)[3]),rep(2,dim(ctrls)[2]*dim(ctrls)[3])))
        PVAL_raw[j]<-ktest$p.value
    }
    
    PVAL_fdr<-p.adjust(PVAL_raw,method='fdr')
    return(list( PVAL_fdr))
}  



wilcox_saver<-function(cells,ctrls){
    
    n.x <-dim(cells)[2]
    n.y <-dim(ctrls)[1]
    n.genes<-dim(cells)[1]

    if(length(dim(cells))==3){
        
        library(foreach)
        library(doSNOW)
        library(parallel)
        
        cluster = makeCluster(5, type = "SOCK")
        registerDoSNOW(cluster)
        getDoParWorkers()
        
        iterations <- dim(cells)[3]
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        
       w= foreach(i=1:dim(cells)[3],.options.snow = opts)%dopar%{
         qq<- apply(cbind(cells[,,i],ctrls[,,i]), 1, FUN = function(z)unlist(wilcox.test(z ~ c(rep(1,dim(cells)[2]),rep(2,dim(ctrls)[2])))[1]))
            return(qq)
       }
       close(pb)
       stopCluster(cluster)
       
       
       mean.w <- Reduce("+", w)/length(w) - n.x*n.y/2
         var.w <- n.x*n.y/12*(n.x+n.y+1) + (1+1/10)/9*apply(simplify2array(w), 1, var)
         mean.z <- sapply(1:n.genes, function(i) (mean.w[i]-sign(mean.w[i])*0.5)/sqrt(var.w[i]))
         mean.p <- sapply(mean.z, function(x) 2*min(pnorm(x), pnorm(x, lower.tail = FALSE)))
         pvaluesrk2 <- cbind(mean.w, mean.p)
         pvaluesrk.adj <- p.adjust(pvaluesrk2[, 2], method = "BH")
        
    } else if(length(dim(cells))==2){
        
        qq<- apply(cbind(cells,ctrls), 1, FUN = function(z)
            wilcox.test(z ~ c(rep(1,dim(cells)[2]),rep(2,dim(ctrls)[2])))$p.value)
        
        pvaluesrk.adj<-p.adjust(qq, method = "BH")
        
        
    }
    
    return(list(pvaluesrk.adj=pvaluesrk.adj))
}    


# x.samp <- t(sapply(1:200, function(i) rgamma(100, 2, 3)))
# 
# w<- apply(x.samp, 1,FUN = function(z) unlist(wilcox.test(z ~ c(rep(1,50),rep(2,50)))[1]))
# 
# unlist(wilcox.test(x.samp[1,] ~ c(rep(1,50),rep(2,50))))[1]

