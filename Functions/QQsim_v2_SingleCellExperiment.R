##############Estimating###########
newQQParams <- function(update = NULL,...) {

  params <- new("QQParams")
  checkmate::assertList(update, null.ok = TRUE)
  update <- c(update, list(...))

  if (length(update) > 0) {
    for (name in names(update)) {
      value <- update[[name]]
      slot( params ,name) <- value
    }
  }
  validObject( params)
  return(  params)
}



QQinitiate<-function(Est_params){
  nGenes <- Est_params@nGenes
  nCells<-ifelse(Est_params@nGroups==1,Est_params@nCells,sum(Est_params@groupCells))
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))


  cells <-  data.frame(Cell = cell.names)
  rownames(cells) <- cell.names
  features <- data.frame(Gene = gene.names)
  rownames(features) <- gene.names

  sim <- SingleCellExperiment(rowData = features, colData = cells, metadata = list(params = params))
  return(sim)
}



setClass("Params",
         contains = "VIRTUAL",
         slots = c(nGenes = "numeric",
                   nCells = "numeric",
                   seed = "numeric"),
         prototype = prototype(nGenes = 10000, nCells = 100,
                               seed = sample(1:1e6, 1)))

setClass("QQParams",
         contains = "Params",
         slots = c(nGroups = "numeric",
                   groupCells = "numeric",
                   mean.shape = "numeric",
                   mean.rate = "numeric",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   beta.loc = "numeric",
                   beta.scale = "numeric",
                   out.prob = "numeric",
                   out.neg.prob= "numeric",
                   out.facLoc = "numeric",
                   out.facScale = "numeric",
                   de.prob = "numeric",
                   de.downProb = "numeric",
                   de.facLoc = "numeric",
                   de.facScale = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   counts.origin="matrix",
                   counts.norm.TC="matrix",
                   Genes.keep = "logical",
                   Beta="numeric",
                   MeanBeta='numeric'),
         prototype = prototype(nGroups = 1,
                               groupCells = 100,
                               mean.rate = 2.109,
                               mean.shape = 1.895,
                               lib.loc = 10.22,
                               lib.scale = 0.3779,
                               beta.loc = 10.22,
                               beta.scale = 0.3779,
                               out.prob = 0.05,
                               out.neg.prob=0,
                               out.facLoc = 2.3184,
                               out.facScale = 0.7462,
                               de.prob = 0.1,
                               de.downProb = 0.5,
                               de.facLoc = 0.1,
                               de.facScale = 0.4,
                               bcv.common = 0.1268,
                               bcv.df = 60,
                               MeanBeta=0.07
         ))
#Param<-new("QQParams")
QQ_Est <- function(counts,Params) {
  Params<-new("QQParams")
  # Normalise for library size and remove all zero genes
  lib.sizes <- colSums(counts)
  lib.med <- median(lib.sizes)
  norm.counts <- t(t(counts) / lib.sizes * lib.med)
  norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

  keep_logit<-rowSums(norm.counts > 0) > 1

  Params@counts.origin<-counts
  Params@counts.norm.TC<-norm.counts
  Params@nGenes<-nrow(counts)
  Params@nCells<-ncol(counts)
  Params@Genes.keep<-keep_logit

  return(Params)
}


winsorize <- function(x, q) {

  checkmate::check_numeric(x, any.missing = FALSE)
  checkmate::check_number(q, lower = 0, upper = 1)

  lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)

  if (diff(lohi) < 0) { lohi <- rev(lohi) }

  x[!is.na(x) & x < lohi[1]] <- lohi[1]
  x[!is.na(x) & x > lohi[2]] <- lohi[2]

  return(x)
}



QQEstLib <- function(Est_params) {
  lib.sizes <- colSums(Est_params@counts.origin)
  fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")

  Est_params@lib.loc<-unname(fit$estimate["meanlog"])
  Est_params@lib.scale<-unname(fit$estimate["sdlog"])

  return(Est_params)
}

QQEst_BETA_MEAN <- function(Est_params,inputBeta=NULL,inputMeanBeta=0.07,inputtrim=0.01) {

  counts<-Est_params@counts.origin


  if(is.null(inputBeta)){
    Trim_temp<-apply((counts),2,function(x){mean(x,trim=inputtrim)})
    Trim_beta<-Trim_temp/mean(Trim_temp)
  }
  else{
    Trim_beta<-inputBeta
  }
  #Est BETA para
  beta_fit <- fitdistrplus::fitdist(Trim_beta/inputMeanBeta, "lnorm")
  #end of BETA est

  norm.counts<-t(t(counts)/Trim_beta/inputMeanBeta)
  means <- rowMeans(norm.counts)
  means <- means[means != 0]
  means <- winsorize(means, q = 0.1)
  fit <- fitdistrplus::fitdist(means, "gamma", method = "mge",gof = "CvM")
  fit$estimate["shape"]
  fit$estimate["rate"]

  if (fit$convergence > 0) {
    warning("Goodness of fit failed, using Method of Moments")
    fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
  }


  Est_params@beta.loc<-unname(beta_fit$estimate["meanlog"])
  Est_params@beta.scale<-unname(beta_fit$estimate["sdlog"])

  Est_params@mean.shape<-unname(fit$estimate["shape"])
  Est_params@mean.rate<-unname(fit$estimate["rate"])

  Est_params@Beta<-Trim_beta
  Est_params@MeanBeta<-inputMeanBeta

  return(Est_params)
}


QQEstOutlier <- function(Est_params) {

  norm.counts<-Est_params@counts.norm.TC

  means <- rowMeans(norm.counts)
  lmeans <- log(means)
  med <- median(lmeans)
  mad <- mad(lmeans)
  bound <- med + 2 * mad
  outs <- which(lmeans > bound)
  prob <- length(outs) / nrow(norm.counts)

  if (length(outs) > 1) {
    facs <- means[outs] / median(means)
    fit <- fitdistrplus::fitdist(facs, "lnorm")}

  Est_params@out.prob<-prob
  Est_params@out.facLoc<-unname(fit$estimate["meanlog"])
  Est_params@out.facScale<- unname(fit$estimate["sdlog"])


  return( Est_params)
}




QQEstBCV <- function(Est_params) {
  # Add dummy design matrix to avoid print statement
  counts<-Est_params@counts.origin
  design <- matrix(1, ncol(counts), 1)
  disps <- edgeR::estimateDisp(counts, design = design)
  Est_params@bcv.common<-0.1 + 0.25 * disps$common.dispersion
  Est_params@bcv.df<-disps$prior.df

  return(Est_params)
}









########Simulating##########
library(bayNorm)
setClass("Sim_params",
         contains = "Params",
         slots = c(
           BaseGeneMean= "numeric",
           OutlierFactor="numeric",
           GeneMean="matrix",
           ExpLibSize='numeric',
           Beta='numeric',
           BCV="matrix",
           CellMeans="matrix",
           TrueCounts="matrix",
           TrueControls="matrix",
           downsample.counts="matrix",
           downsample.controls="matrix",
           BaseCellMeans='matrix',
           Group_DEfacs_list='list',
           Group_means_gene_list='list',
           GroupInd='numeric',
           BaseCellMeans_group= "matrix"
         ))


getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale,Prior=NULL) {


  if(is.null(Prior)){
    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
  }


  else{

    probb<-ifelse(Prior<  quantile(Prior,0.05),0,sel.prob)

    is.selected <- as.logical(rbinom(n.facs, 1, probb))
  }


  n.selected <- sum(is.selected)
  dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
  facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)

  # Reverse directions for factors that are less than one
  dir.selected[facs.selected < 1 & dir.selected == -1] <- 1
  dir.selected[facs.selected < 1 & dir.selected == 1] <- -1
  factors <- rep(1, n.facs)
  factors[is.selected] <- facs.selected ^ dir.selected

  return(factors)
}


QQSimGeneMeans <- function(SCE,Est_params,RandomPCR=F) {

  cell.names <- colData(SCE)
  gene.names <- rowData(SCE)

  nGenes <- Est_params@nGenes

  nCells<-ifelse(Est_params@nGroups==1,Est_params@nCells,sum(Est_params@groupCells))
  #nCells<-Est_params@nCells

  mean.shape <- Est_params@mean.shape
  mean.rate <-Est_params@mean.rate
  out.prob <- Est_params@out.prob
  out.facLoc <- Est_params@out.facLoc
  out.facScale <-Est_params@out.facScale
  out.neg.prob<-Est_params@out.neg.prob
  # Simulate base gene means
  base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

  if(RandomPCR){
    library(foreach)
    Outlier_mat<-matrix(nrow=nGenes,ncol=nCells)
    means.gene<-foreach(i=1:nCells,.combine=cbind)%do%{
      #Randomize outlier.facs for different cells:
      if(out.prob==0){outlier.facs<-rep(1,nGenes)}else{
        outlier.facs <- getLNormFactors(nGenes, out.prob,out.neg.prob, out.facLoc,out.facScale,Prior=NULL)

      }

      Outlier_mat[,i]<-outlier.facs

      # Add expression outliers
      median.means.gene <- median(base.means.gene)
      outlier.means <- median.means.gene * Outlier_mat[,i]
      is.outlier <- outlier.facs != 1
      means.gene <- base.means.gene
      means.gene[is.outlier] <- outlier.means[is.outlier]
      return(means.gene)
    }

    colnames(Outlier_mat) <- cell.names
    rownames(Outlier_mat) <- gene.names
    assays(SCE)$Outlier_mat <-  Outlier_mat

  } else{
    if(out.prob==0){
      outlier.facs<-rep(1,nGenes)
    }else{
      outlier.facs <- getLNormFactors(nGenes, out.prob, out.neg.prob, out.facLoc,out.facScale,Prior=NULL)

    }

    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    means.gene<-matrix(means.gene,ncol=nCells ,nrow=nGenes,byrow=F)

    rowData(SCE)$OutlierFactor <- outlier.facs
  }

  colnames(means.gene) <- cell.names$Cell
  rownames(means.gene) <- gene.names$Gene

  rowData(SCE)$BaseGeneMean <- base.means.gene
  assays(SCE)$Single_Genemean_mat <- means.gene
  return(SCE)
}



QQSimGroupDE <- function(SCE,Est_params) {

  nGenes <- Est_params@nGenes
  nGroups <- Est_params@nGroups
  de.prob <- Est_params@de.prob
  de.downProb <- Est_params@de.downProb
  de.facLoc <- Est_params@de.facLoc
  de.facScale <- Est_params@de.facScale
  groupCells<-Est_params@groupCells

  BaseGeneMean.out<- rowMeans(assays(SCE)$Single_Genemean_mat)

  rowData(SCE)$GroupInd<-NULL
  print('de.prob, de.downProb, de.facLoc and de.facScale must be vectors')

  for (idx in seq_len(nGroups)) {
    de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx],de.facLoc[idx], de.facScale[idx])
    group.means.gene <-BaseGeneMean.out * de.facs

    rowData(SCE)[[paste0("GeneMeanGroup", idx)]] <- group.means.gene
    rowData(SCE)[[paste0("DEFacGroup", idx)]]<-de.facs

  }

  colData(SCE)$GroupInd<-rep(seq(1,nGroups),groupCells)
  #colData(SCE)$Group<-paste0("Group", colData(SCE)$GroupInd)
  colData(SCE)$Group<-colData(SCE)$GroupInd

  return(SCE)
}




#qq<-Sim_List_Input$Sim_params@GeneMean
QQSimBETA <- function(SCE, Est_params) {
  nGroups <- Est_params@nGroups
  GroupInd<-colData(SCE)$GroupInd
  groups <- colData(SCE)$Group
  group.names <- unique(groups)

  Est_params@nCells<-ifelse(Est_params@nGroups==1,Est_params@nCells,sum(Est_params@groupCells))
  nCells<-Est_params@nCells

  lib.loc <- Est_params@lib.loc
  lib.scale <- Est_params@lib.scale
  beta.loc <- Est_params@beta.loc
  beta.scale <- Est_params@beta.scale


  counts<-Est_params@counts.origin


  exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
  colData(SCE)$ExpLibSize <- exp.lib.sizes

  ##Beta: for simulation without real data###
  if(length(Est_params@Beta)==0 || is.null(Est_params@Beta)){
    if(nGroups==1){
      temp<-rlnorm(nCells,meanlog=beta.loc,sdlog=beta.scale)
      colData(SCE)$Beta<-temp/mean(temp)
    }else{

      library(foreach)
      colData(SCE)$Beta <-foreach(idx= 1:nGroups,.combine=c)%do%{

        temp<-rlnorm(length(which(GroupInd==idx)), meanlog=beta.loc,sdlog=beta.scale)
        return(temp/mean(temp))

      }

    }

  } else{
    colData(SCE)$Beta<-Est_params@Beta
  }
  return(SCE)
}




QQSimGroupCellMeans <- function(SCE, Est_params,Effect=F) {
  cell.names <- colData(SCE)$Cell
  gene.names <- rowData(SCE)$Gene

  nGroups <- Est_params@nGroups
  exp.lib.sizes <-colData(SCE)$ExpLibSize
  GroupInd<-colData(SCE)$GroupInd
  groups <- colData(SCE)$Group
  group.names <- unique(groups)
  group.means.gene<- rowData(SCE)[, paste0("GeneMeanGroup", group.names)]

  Group_Genemean_mat<-as.matrix(group.means.gene[, groups])


  colnames(Group_Genemean_mat) <- cell.names
  rownames(Group_Genemean_mat) <- gene.names
  assays(SCE)$Group_Genemean_mat <- Group_Genemean_mat

  return(SCE)
}



QQSimBCVMeans <- function(SCE,Est_params,SIZE_MU_Dependence=T,shape=1,scale=1) {


  cell.names <- colData(SCE)$Cell
  gene.names <- rowData(SCE)$Gene
  nGenes <- Est_params@nGenes
  #nCells <- Est_params@nCells
  nCells<-ifelse(Est_params@nGroups==1,Est_params@nCells,sum(Est_params@groupCells))
  bcv.common <- Est_params@bcv.common
  bcv.df <- Est_params@bcv.df

  if(Est_params@nGroups==1){
    GeneMean_mat<-assays(SCE)$Single_Genemean_mat


  }else{
    GeneMean_mat <-assays(SCE)$Group_Genemean_mat
  }

  if(SIZE_MU_Dependence){
    if (is.finite(bcv.df)) {
      bcv <- (bcv.common + (1 / sqrt(GeneMean_mat))) *
        sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    } else {
      warning("'bcv.df' is infinite. This parameter will be ignored.")
      bcv <- (bcv.common + (1 / sqrt(GeneMean_mat)))
    }

  } else{

    phi=rgamma(nGenes,shape=shape,scale=scale)
    bcv<-1/sqrt(phi)
    bcv<-matrix(bcv,ncol=nCells ,nrow=nGenes,byrow=F)
  }

  rownames(bcv)<-gene.names
  colnames(bcv)<-cell.names
  assays(SCE)$BCV<- bcv


  Gamma_Means_mat <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2), scale = GeneMean_mat * (bcv ^ 2)),nrow = nGenes, ncol = nCells)

  true.counts <- matrix(rpois(nGenes * nCells, lambda =  Gamma_Means_mat ),nrow = nGenes, ncol = nCells)
  colnames(true.counts) <- cell.names
  rownames(true.counts) <- gene.names
  assays(SCE)$TrueCounts<- true.counts
  return(SCE)
}

library(foreach)
QQSimBinomial<-function(SCE,Est_params){
  nGroups <- Est_params@nGroups
  GroupInd<-colData(SCE)$GroupInd
  groups <- colData(SCE)$Group
  group.names <- unique(groups)

  cell.names <- colData(SCE)$Cell
  gene.names <- rowData(SCE)$Gene

  nGenes <- Est_params@nGenes

  nCells<-ifelse(Est_params@nGroups==1,Est_params@nCells,sum(Est_params@groupCells))
  Beta<-colData(SCE)$Beta

  if(Est_params@nGroups==1){
      bb<-colData(SCE)$Beta*Est_params@MeanBeta
      bb[bb>=1]<-max(bb[bb<1])

    downsample.counts <-bayNorm::DownSampling(assays(SCE)$TrueCounts ,bb)
  }else{

    downsample.counts <-foreach(idx= 1:nGroups,.combine=cbind)%do%{
        bb<-colData(SCE)$Beta[which(GroupInd==idx)]*Est_params@MeanBeta[idx]
        bb[bb>=1]<-max(bb[bb<1])

      temp<-assays(SCE)$TrueCounts [,which(GroupInd==idx)]
      temp.down<- bayNorm::DownSampling(temp ,BETA_vec =bb)
      return(temp.down)
    }

  }

  colnames( downsample.counts) <- cell.names
  rownames( downsample.counts) <- gene.names
  assays(SCE)$counts<-downsample.counts
  return(SCE)
}


trytt2<-function(counts,inputBeta=NULL,inputMeanBeta=0.07,inputtrim=0.01){

  ##First step##:
  Est_params<-QQ_Est(as.matrix(counts))
  Est_params<-QQEstLib(Est_params)
  Est_params<-QQEst_BETA_MEAN(Est_params,inputBeta=inputBeta,inputMeanBeta=inputMeanBeta,inputtrim=inputtrim)
  Est_params<-QQEstOutlier(Est_params)
  system.time(Est_params<-QQEstBCV(Est_params))

  #Simulation
  SCE<-QQinitiate(Est_params)
  SCE<-QQSimGeneMeans(SCE,Est_params)
  SCE<-QQSimBETA(SCE, Est_params)
  SCE<-QQSimBCVMeans(SCE,Est_params)
  #SCE<-QQSimTrueCounts(SCE, Est_params)
  SCE<-QQSimBinomial(SCE,Est_params)

  return(list(SCE=SCE,Est_params=Est_params))
}


# #debug##########
# data("sc_example_counts")
# params <- splatEstimate(sc_example_counts)
#
# ##First step##:
# Est_params<-QQ_Est(as.matrix(sc_example_counts))
# Est_params<-QQEstLib(Est_params)
# Est_params<-QQEst_BETA_MEAN(Est_params,inputBeta=NULL,inputMeanBeta=0.05,inputtrim=0.05)
# Est_params<-QQEstOutlier(Est_params)
# system.time(Est_params<-QQEstBCV(Est_params))
#
# #Simulation
# SCE<-QQinitiate(Est_params)
# SCE<-QQSimGeneMeans(SCE,Est_params)
# SCE<-QQSimBETA(SCE, Est_params)
# SCE<-QQSimBCVMeans(SCE,Est_params)
# #SCE<-QQSimTrueCounts(SCE, Est_params)
# SCE<-QQSimBinomial(SCE,Est_params)
