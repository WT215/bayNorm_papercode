######################### load functions##########
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



addFeatureStats <- function(sce, value = c("counts", "cpm", "tpm", "fpkm"),
                            log = FALSE, offset = 1, no.zeros = FALSE) {

    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertLogical(log)
    checkmate::assertNumber(offset, lower = 0)
    checkmate::assertLogical(no.zeros)
    value <- match.arg(value)

    switch(value,
           counts = {
               values = BiocGenerics::counts(sce)
               suffix <- "Counts"
           },
           cpm = {
               values = SingleCellExperiment::cpm(sce)
               suffix <- "CPM"
           },
           tpm = {
               values = SingleCellExperiment::tpm(sce)
               suffix <- "TPM"
           },
           fpkm = {
               values = SummarizedExperiment::assays(sce)$fpkm
               suffix <- "FPKM"
           }
    )

    if (no.zeros) {
        values[values == 0] <- NA
        suffix = paste0(suffix, "No0")
    }

    if (log) {
        values = log2(values + offset)
        suffix = paste0("Log", suffix)
    }

    mean.str <- paste0("Mean", suffix)
    var.str  <- paste0("Var",  suffix)
    cv.str   <- paste0("CV",   suffix)
    med.str  <- paste0("Med",  suffix)
    mad.str  <- paste0("MAD",  suffix)

    rowData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    rowData(sce)[, var.str]  <- matrixStats::rowVars(values, na.rm = TRUE)
    rowData(sce)[, cv.str]   <- sqrt(rowData(sce)[, var.str]) /
        rowData(sce)[, mean.str]
    rowData(sce)[, med.str]  <- matrixStats::rowMedians(values, na.rm = TRUE)
    rowData(sce)[, mad.str]  <- matrixStats::rowMads(values, na.rm = TRUE)

    return(sce)
}
rbindMatched <- function(df1, df2) {
    common.names <- intersect(colnames(df1), colnames(df2))
    combined <- rbind(df1[, common.names], df2[, common.names])

    return(combined)
}

library(grid)
#####plot_DROPOUT########
cbPaletteee <- c("#999999", "#E69F00", "#56B4E9")
#cbPaletteee<-c("#F8766D" ,"#00BFC4","#619CFF")
library(ggplot2)
plot_DROPOUT<-function(listused_N1,MAIN='',CAPTION='',legendpointsize=1,legend_key_size=1,subtitle=''){

library(foreach)
DROPOUT_DAT<-foreach(i=1:length(listused_N1),.combine=rbind)%do%{

    if(names(listused_N1)[i]=='Real') {colll<-cbPaletteee[1]}else if(names(listused_N1)[i]=='Binomial'){
        colll<-cbPaletteee[2]
    } else if(names(listused_N1)[i]=='Scaled raw'){
        colll<-cbPaletteee[3]
    } else{colll<-cbPaletteee[1]}

  dropout<-apply(listused_N1[[i]],1,function(x){length(which(x==0))/length(x)})
  meann<-rowMeans(listused_N1[[i]])
  qqq<-cbind(dropout,meann,rep(names(listused_N1)[i],length(dropout)),rep(colll,length(dropout)))
  return(qqq)
}
DROPOUT_DAT<-as.data.frame(DROPOUT_DAT)
colnames(DROPOUT_DAT)<-c('Dropout rate','Mean expression','Dataset','Colour')
DROPOUT_DAT$`Dropout rate`<-as.numeric(as.character(DROPOUT_DAT$`Dropout rate`))
DROPOUT_DAT$`Mean expression`<-as.numeric(as.character(DROPOUT_DAT$`Mean expression`))
DROPOUT_DAT$Dataset<-factor(DROPOUT_DAT$Dataset,levels=unique(DROPOUT_DAT$Dataset))
DROPOUT_DAT$Colour<-factor(DROPOUT_DAT$Colour,levels=unique(DROPOUT_DAT$Colour))



sces<-listused_N1
colours <- scales::hue_pal()(length(sces))

theoline_ref_name<-names(listused_N1)[1]

theoline<-data.frame(xx =sort(DROPOUT_DAT$`Mean expression`[which(DROPOUT_DAT$Dataset==theoline_ref_name)]), yy=exp(-sort(DROPOUT_DAT$`Mean expression`[which(DROPOUT_DAT$Dataset==theoline_ref_name)])))

mean.zeros <- ggplot() +
    geom_line(data=theoline,aes(x=xx,y=yy, lty = 'exp(-mean expression)'))+
    geom_point(data=DROPOUT_DAT,aes_string(x = DROPOUT_DAT$`Mean expression`, y = DROPOUT_DAT$`Dropout rate`,colour = DROPOUT_DAT$Colour),size = point.size, alpha = point.alpha,shape=46) +
  scale_x_log10(labels = scales::comma) +
  xlab("Mean expression") +
  ylab("Dropout rates") +
    scale_color_manual(values=as.character(unique(DROPOUT_DAT$Colour)),labels=names(listused_N1))+
    scale_linetype_manual(values=1,'Black line')+
  labs(x = "Mean expression",y="Dropout rates",fill='Dataset',colour='Dataset',caption=CAPTION,title=MAIN,subtitle=subtitle)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = 10),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize),legend.key.size = unit(legend_key_size,"line") ,legend.position = c(0.8,0.9))+
    guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))
return(mean.zeros)
}
#

#


#####mean var#####

plot_MEANVAR<-function(sces,MAIN='',CAPTION='',legend.position='top'){

    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- scater::calculateCPM(sce, use_size_factors = FALSE)
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$total_features
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        sces[[name]] <- sce
    }

    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }
    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)


    #colours <- scales::hue_pal()(length(sces))
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbbPalette<-c("black", "red", "red2", "red4", "green", "green2", "green4", "blue", "blue2", "blue4")
    colours <- cbbPalette[seq(1,length(sces))]



    z.gene <- ggplot(features,
                     aes_string(x = "Dataset", y = "pct_dropout_counts",
                                colour = "Dataset")) +
        geom_violin(show.legend=F)+
        #geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage dropout rates per gene") +
        ggtitle("Distribution of dropout rates per gene") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank())

    z.cell <- ggplot(cells,
                     aes_string(x = "Dataset", y = "PctZero",
                                colour = "Dataset")) +
        geom_violin(show.legend=F)+
        #geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage dropout rates per cell") +
        ggtitle("Distribution of dropout rates per cell") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank())


    library(foreach)
    Hicks_line1<- foreach(i=1:length(sces),.combine=rbind)%do%{
        xx<-features$MeanLogCPM[features$Dataset==names(sces)[i]]
        yy<-features$VarLogCPM[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =25)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line1)[5]<-'Dataset'
    
    
    mean.var <- #ggplot()+
        ggplot(features,aes_string(x = "MeanLogCPM", y = "VarLogCPM",colour = "Dataset", fill = "Dataset")) +
        #geom_point(size = point.size, alpha = point.alpha,shape=46) +
    geom_line(data=Hicks_line1,aes_string(x=Hicks_line1$x,y=Hicks_line1$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8,show.legend = T)+
        #geom_smooth(size = linewidth, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))+
        # xlab(expression(paste("Mean", log[2], "(CPM + 1)"))) +
        # ylab(expression(paste("Variance", log[2], "(CPM + 1)"))) +
        xlab("Mean expression") +
        ylab("Variance of gene expression") +
        labs(caption=CAPTION,colour='',fill='')+
        ggtitle("Mean-variance relationship",subtitle = MAIN) +
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(legend_key_size,"line"),legend.position =legend.position,legend.background=element_blank(),legend.key =element_blank())

    mean.var_v2 <- ggplotGrob(mean.var)
    # Check out the grobs
    library(grid)
    names.grobs <- grid.ls(grid.force( mean.var_v2), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        mean.var_v2 <- editGrob(grid.force( mean.var_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =colours [i]))
    }
    
    #as_ggplot(mean.var_v2)
    
  ###zeros trend  
    theoline_ref_name<-names(sces)[1]
    theoline<-data.frame(xx =features$mean_counts[which(features$Dataset== theoline_ref_name)], yy=exp(-features$mean_counts[which(features$Dataset== theoline_ref_name)]))
    
    
    library(foreach)
    Hicks_line2<- foreach(i=1:length(sces),.combine=rbind)%do%{
        yy<-features$pct_dropout_by_counts[features$Dataset==names(sces)[i]]/100
        xx<-features$mean_counts[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =100)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line2)[5]<-'Dataset'

    mean.zeros <- ggplot() +

        # geom_point(data=features,
        #            aes_string(x = "MeanCounts",
        #                       y = features$pct_dropout_by_counts/100,colour = "Dataset",
        #                       fill = "Dataset"),
        #            size = point.size,
        #            alpha = point.alpha,shape=46) +
        
        # geom_smooth(data=features,
        #             aes_string(x = "MeanCounts",
        #                        y = features$pct_dropout_by_counts/100,colour = "Dataset",fill = "Dataset"),
        #             size = linewidth,
        #             alpha = point.alpha) +
    


        geom_line(data=Hicks_line2,aes_string(x=Hicks_line2$x,y=Hicks_line2$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8,show.legend=F)+
        geom_line(aes(x=theoline$xx,y=theoline$yy, linetype = 'exp(-mean expression)'),size = linewidth.exp,alpha = point.alpha)+
        scale_x_log10(labels = scales::comma) +
        scale_colour_manual(values = colours) +
        scale_linetype_manual(values='dashed','Dashed line')+
        #scale_fill_manual(values = colours) +
        guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))+
        labs(x = "Mean expression",y="Dropout rates",fill='',colour='',caption=CAPTION,title="Mean-dropout rates relationship",subtitle='')+
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
            theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(legend_key_size,"line"),legend.position = legend.position,legend.background=element_blank(),legend.key =element_blank())

    
    
    # mean.zeros_v2 <- ggplotGrob(mean.zeros)
    # # Check out the grobs
    # library(grid)
    # names.grobs <- grid.ls(grid.force( mean.zeros_v2), print=F)$name 
    # labels <- names.grobs[which(grepl("^label", names.grobs))]
    # 
    # for(i in seq_along(labels)[-1]) {
    #     mean.zeros_v2 <- editGrob(grid.force(mean.zeros_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =c("#FFFFFF66",colours) [i]))
    # }
    # #as_ggplot( mean.zeros_v2)
    # g <- ggplot_build(mean.zeros)
    # unique(g$data[[2]]$colour)

    return(list(mean.var=mean.var_v2,z.gene=z.gene,z.cell=z.cell,mean.zeros=mean.zeros))
}
##########difff##############
plot_MEANVAR_diff<-function(sces,MAIN='',CAPTION=''){

    ref<-1
    ref.dim <- dim(sces[[ref]])

    for (name in names(sces)) {
        sce <- sces[[name]]
        if (!identical(dim(sce), ref.dim)) {
            stop("all datasets in 'sces' must have the same dimensions")
        }
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- scater::calculateCPM(sce, use_size_factors  = FALSE)
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$total_features
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        rowData(sce)$RankCounts <- rank(rowData(sce)$mean_counts)
        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]

    ref.means <- sort(rowData(ref.sce)$MeanLogCPM)
    ref.vars <- sort(rowData(ref.sce)$VarLogCPM)
    ref.libs <- sort(colData(ref.sce)$total_counts)
    ref.z.gene <- sort(rowData(ref.sce)$pct_dropout_counts)
    ref.z.cell <- sort(colData(ref.sce)$PctZero)

    ref.rank.ord <- order(rowData(ref.sce)$RankCounts)
    ref.vars.rank <- rowData(ref.sce)$VarLogCPM[ref.rank.ord]
    ref.z.gene.rank <- rowData(ref.sce)$pct_dropout_counts[ref.rank.ord]

    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$RefRankMeanLogCPM <- ref.means[
            rank(rowData(sce)$MeanLogCPM)]
        rowData(sce)$RankDiffMeanLogCPM <- rowData(sce)$MeanLogCPM -
            rowData(sce)$RefRankMeanLogCPM
        rowData(sce)$RefRankVarLogCPM <- ref.vars[rank(rowData(sce)$VarLogCPM)]
        rowData(sce)$RankDiffVarLogCPM <- rowData(sce)$VarLogCPM -
            rowData(sce)$RefRankVarLogCPM
        colData(sce)$RefRankLibSize <- ref.libs[rank(colData(sce)$total_counts)]
        colData(sce)$RankDiffLibSize <- colData(sce)$total_counts -
            colData(sce)$RefRankLibSize
        rowData(sce)$RefRankZeros <- ref.z.gene[rank(
            rowData(sce)$pct_dropout_counts)]
        rowData(sce)$RankDiffZeros <- rowData(sce)$pct_dropout_counts -
            rowData(sce)$RefRankZeros
        colData(sce)$RefRankZeros <- ref.z.cell[rank(
            colData(sce)$PctZero)]
        colData(sce)$RankDiffZeros <- colData(sce)$PctZero -
            colData(sce)$RefRankZeros

        rowData(sce)$MeanRankVarDiff <- rowData(sce)$VarLogCPM -
            ref.vars.rank[rowData(sce)$RankCounts]
        rowData(sce)$MeanRankZerosDiff <- rowData(sce)$pct_dropout_counts -ref.z.gene.rank[rowData(sce)$RankCounts]

        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]
    sces[[ref]] <- NULL

    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])

    
    #colours <- scales::hue_pal()(length(sces)+1)
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbbPalette<-c("black", "red", "red2", "red4", "green", "green2", "green4", "blue", "blue2", "blue4")
    colours <- cbbPalette[seq(1,length(sces)+1)]

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }

    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)


    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }

    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)



    #plot(features$RankCounts[features$Dataset=='Binomial'],features$MeanRankVarDiff[features$Dataset=='Binomial'])
    
    
    library(foreach)
    Hicks_line_meanvar<- foreach(i=1:length(sces),.combine=rbind)%do%{
        yy<-features$MeanRankVarDiff[features$Dataset==names(sces)[i]]
        xx<-features$RankCounts[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =25)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line_meanvar)[5]<-'Dataset'
    

    mean.var <- ggplot(features,
                       aes_string(x = "RankCounts", y = "MeanRankVarDiff",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,size=hline_size,linetype='dotted') +
        #geom_point(size = point.size, alpha = point.alpha,shape=46) +
        geom_smooth(size = linewidth,show.legend=F)+
        #geom_line(data=Hicks_line_meanvar,aes_string(x=Hicks_line_meanvar$x,y=Hicks_line_meanvar$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8)+
        scale_colour_manual(values = colours[-1]) +
        scale_fill_manual(values = colours[-1]) +
        xlab("Expression rank") +
        ylab(expression(paste("Difference in variance "))) +
        #ylab(expression(paste("Difference in variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Difference in mean-variance relationship") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom')



    z.gene <- ggplot(features,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,linewidth=hline_size,linetype='dotted') +
        geom_boxplot(show.legend=F,outlier.shape=46) +
        scale_colour_manual(values = colours[-1]) +
        ylab(paste("Difference in dropout rates per gene")) +
        ggtitle("Difference in dropout rates per gene") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank())

    z.cell <- ggplot(cells,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,linewidth=hline_size,linetype='dotted') +
        geom_boxplot(show.legend=F,outlier.shape=46) +
        scale_colour_manual(values = colours[-1]) +
        ylab(paste("Difference in dropout rates per cell")) +
        ggtitle("Difference in dropout rates per cell") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank())

    library(foreach)
    Hicks_line_meanzeros<- foreach(i=1:length(sces),.combine=rbind)%do%{
        yy<-features$MeanRankZerosDiff[features$Dataset==names(sces)[i]]
        xx<-features$RankCounts[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =25)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line_meanzeros)[5]<-'Dataset'
    
    
    # yy<-features$MeanRankZerosDiff[features$Dataset==names(sces)[2]]
    # xx<-features$RankCounts[features$Dataset==names(sces)[2]]
    # plot(xx,yy)
    
    
    mean.zeros <- ggplot(features,aes_string(x = "RankCounts", y = "MeanRankZerosDiff", colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,size=hline_size,linetype='dotted') +
        #geom_point(size = point.size, alpha = point.alpha,shape=46) +
        geom_smooth(size = linewidth,show.legend=F)+
        #geom_line(data=Hicks_line_meanzeros,aes_string(x=Hicks_line_meanzeros$x,y=Hicks_line_meanzeros$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8)+
        scale_colour_manual(values = colours[-1]) +
        scale_fill_manual(values = colours[-1]) +
        xlab("Expression rank") +
        ylab("Difference in dropout rates") +
        ggtitle("Difference in mean-dropout rates relationship") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom')

    return(list(mean.var=mean.var,z.gene=z.gene,z.cell=z.cell,mean.zeros=mean.zeros))
    #return(list(z.gene=z.gene,z.cell=z.cell))
}#end



####meanvar v2 for main figure##########
library(gridExtra)
library(ggpubr)

plot_MEANVAR_v2<-function(sces,MAIN='',CAPTION='',legend.position='top'){
    
    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- scater::calculateCPM(sce, use_size_factors = FALSE)
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$total_features
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        sces[[name]] <- sce
    }
    
    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])
    
    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }
    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)
    
    

    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbbPalette<-c("black", "red", "red2", "red4", "green", "green2", "green4", "blue", "blue2", "blue4")
    colours <- cbbPalette[seq(1,length(sces))]
    
    z.gene <- ggplot(features,
                     aes_string(x = "Dataset", y = "pct_dropout_counts",
                                colour = "Dataset")) +
        geom_violin(show.legend=F)+
        #geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage dropout rates per gene") +
        ggtitle("Distribution of dropout rates per gene") +
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank(),legend.background=element_blank(),legend.key =element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank())
    
    
    z.gene_v2 <- ggplotGrob(z.gene)
    # Check out the grobs
    library(grid)
    names.grobs <- grid.ls(grid.force(z.gene_v2), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        z.gene_v2 <- editGrob(grid.force(z.gene_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =colours [i]))
    }
    #z.gene_v2<-as_ggplot(z.gene_v2)
    
    z.cell <- ggplot(cells,
                     aes_string(x = "Dataset", y = "PctZero",
                                colour = "Dataset")) +
        geom_violin(show.legend=T)+
        #geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage dropout rates per cell") +
        ggtitle("Distribution of dropout rates per cell") +
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        labs(colour='')+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='none',axis.title.x=element_blank(),legend.background=element_blank(),legend.key =element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(),legend.direction = "horizontal")
    
    z.cell_v2 <- ggplotGrob(z.cell)
    # Check out the grobs
    library(grid)
    names.grobs <- grid.ls(grid.force( z.cell_v2), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        z.cell_v2 <- editGrob(grid.force( z.cell_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =colours [i]))
    }
    #z.cell_v2<-as_ggplot( z.cell_v2)
    
    
    library(foreach)
    mean.var <- #ggplot()+
        ggplot(features,aes_string(x = "MeanLogCPM", y = "VarLogCPM",colour = "Dataset", fill = "Dataset")) +
        geom_point(size = point.size, alpha = point.alpha,shape=46) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))+
        # xlab(expression(paste("Mean", log[2], "(CPM + 1)"))) +
        # ylab(expression(paste("Variance", log[2], "(CPM + 1)"))) +
        xlab("Mean expression") +
        ylab("Variance of gene expression") +
        labs(caption=CAPTION,fill='',colour='')+
        ggtitle("Mean-variance relationship") +
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(legend_key_size,"line"),legend.position =legend.position,legend.background=element_blank(),legend.key =element_blank())
    
    mean.var_v2 <- ggplotGrob(mean.var)
    # Check out the grobs
    library(grid)
    names.grobs <- grid.ls(grid.force( mean.var_v2), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        mean.var_v2 <- editGrob(grid.force( mean.var_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =colours [i]))
    }
    #mean.var_v2<-as_ggplot( mean.var_v2)
    
    
    theoline_ref_name<-names(sces)[1]
    theoline<-data.frame(xx =features$mean_counts[which(features$Dataset== theoline_ref_name)], yy=exp(-features$mean_counts[which(features$Dataset== theoline_ref_name)]))
    
    
    
    mean.zeros <- ggplot() +
        
        geom_point(data=features,
                   aes_string(x = "MeanCounts",
                              y = features$pct_dropout_by_counts/100,colour = "Dataset",
                              fill = "Dataset"),
                   size = point.size,
                   alpha = point.alpha,shape=46,show.legend=F) +

        geom_line(aes(x=theoline$xx,y=theoline$yy, linetype = 'exp(-mean expression)'),size = linewidth.exp,alpha = point.alpha)+
        scale_x_log10(labels = scales::comma) +
        scale_colour_manual(values = colours) +
        scale_linetype_manual(values='dashed','Dashed line')+
        #scale_fill_manual(values = colours) +
        guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))+
        labs(x = "Mean expression",y="Dropout rates",fill='Dataset',colour='Dataset',caption=CAPTION,title="Mean-dropout rates relationship",subtitle='')+
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(legend_key_size,"line"),legend.position = legend.position,legend.background=element_blank(),legend.key =element_blank())
    
    mean.zeros_v2 <- ggplotGrob(mean.zeros)
    # Check out the grobs
    library(grid)
    names.grobs <- grid.ls(grid.force( mean.zeros_v2), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        mean.zeros_v2 <- editGrob(grid.force( mean.zeros_v2), gPath(labels[i]), grep = TRUE,  gp = gpar(col =colours [i]))
    }
    #mean.zeros_v2<-as_ggplot( mean.zeros_v2)
    

    return(list(mean.var=mean.var_v2,z.gene=z.gene_v2,z.cell=z.cell_v2,mean.zeros=mean.zeros_v2))
}





####diffff v2########
plot_MEANVAR_diff_v2<-function(sces,MAIN='',CAPTION=''){
    
    ref<-1
    ref.dim <- dim(sces[[ref]])
    
    for (name in names(sces)) {
        sce <- sces[[name]]
        if (!identical(dim(sce), ref.dim)) {
            stop("all datasets in 'sces' must have the same dimensions")
        }
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- scater::calculateCPM(sce, use_size_factors  = FALSE)
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$total_features
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        rowData(sce)$RankCounts <- rank(rowData(sce)$mean_counts)
        sces[[name]] <- sce
    }
    
    ref.sce <- sces[[ref]]
    
    ref.means <- sort(rowData(ref.sce)$MeanLogCPM)
    ref.vars <- sort(rowData(ref.sce)$VarLogCPM)
    ref.libs <- sort(colData(ref.sce)$total_counts)
    ref.z.gene <- sort(rowData(ref.sce)$pct_dropout_counts)
    ref.z.cell <- sort(colData(ref.sce)$PctZero)
    
    ref.rank.ord <- order(rowData(ref.sce)$RankCounts)
    ref.vars.rank <- rowData(ref.sce)$VarLogCPM[ref.rank.ord]
    ref.z.gene.rank <- rowData(ref.sce)$pct_dropout_counts[ref.rank.ord]
    
    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$RefRankMeanLogCPM <- ref.means[
            rank(rowData(sce)$MeanLogCPM)]
        rowData(sce)$RankDiffMeanLogCPM <- rowData(sce)$MeanLogCPM -
            rowData(sce)$RefRankMeanLogCPM
        rowData(sce)$RefRankVarLogCPM <- ref.vars[rank(rowData(sce)$VarLogCPM)]
        rowData(sce)$RankDiffVarLogCPM <- rowData(sce)$VarLogCPM -
            rowData(sce)$RefRankVarLogCPM
        colData(sce)$RefRankLibSize <- ref.libs[rank(colData(sce)$total_counts)]
        colData(sce)$RankDiffLibSize <- colData(sce)$total_counts -
            colData(sce)$RefRankLibSize
        rowData(sce)$RefRankZeros <- ref.z.gene[rank(
            rowData(sce)$pct_dropout_counts)]
        rowData(sce)$RankDiffZeros <- rowData(sce)$pct_dropout_counts -
            rowData(sce)$RefRankZeros
        colData(sce)$RefRankZeros <- ref.z.cell[rank(
            colData(sce)$PctZero)]
        colData(sce)$RankDiffZeros <- colData(sce)$PctZero -
            colData(sce)$RefRankZeros
        
        rowData(sce)$MeanRankVarDiff <- rowData(sce)$VarLogCPM -
            ref.vars.rank[rowData(sce)$RankCounts]
        rowData(sce)$MeanRankZerosDiff <- rowData(sce)$pct_dropout_counts -ref.z.gene.rank[rowData(sce)$RankCounts]
        
        sces[[name]] <- sce
    }
    
    ref.sce <- sces[[ref]]
    sces[[ref]] <- NULL
    
    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])
    
    
    #colours <- scales::hue_pal()(length(sces)+1)
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbbPalette<-c("black", "red", "red2", "red4", "green", "green2", "green4", "blue", "blue2", "blue4")
    colours <- cbbPalette[seq(1,length(sces)+1)]
    
    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }
    
    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)
    
    
    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }
    
    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)
    
    
    
    #plot(features$RankCounts[features$Dataset=='Binomial'],features$MeanRankVarDiff[features$Dataset=='Binomial'])
    
    
    library(foreach)
    Hicks_line_meanvar<- foreach(i=1:length(sces),.combine=rbind)%do%{
        yy<-features$MeanRankVarDiff[features$Dataset==names(sces)[i]]
        xx<-features$RankCounts[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =25)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line_meanvar)[5]<-'Dataset'
    
    
    mean.var <- ggplot(features,
                       aes_string(x = "RankCounts", y = "MeanRankVarDiff",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,size=hline_size,linetype='dotted') +
        #geom_point(size = point.size, alpha = point.alpha,shape=46) +
        geom_smooth(size = linewidth)+
        #geom_line(data=Hicks_line_meanvar,aes_string(x=Hicks_line_meanvar$x,y=Hicks_line_meanvar$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8)+
        scale_colour_manual(values = colours[-1]) +
        scale_fill_manual(values = colours[-1]) +
        xlab("Expression rank") +
        ylab(expression(paste("Difference in variance "))) +
        #ylab(expression(paste("Difference in variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Difference in mean-variance relationship") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom')
    
    
    
    z.gene <- ggplot(features,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,linewidth=hline_size,linetype='dotted') +
        geom_boxplot(show.legend=F,outlier.shape=46) +
        scale_colour_manual(values = colours[-1]) +
        ylab(paste("Difference in dropout rates per gene")) +
        ggtitle("Difference in dropout rates per gene") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank())
    
    z.cell <- ggplot(cells,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,linewidth=hline_size,linetype='dotted') +
        geom_boxplot(show.legend=F,outlier.shape=46) +
        scale_colour_manual(values = colours[-1]) +
        ylab(paste("Difference in dropout rates per cell")) +
        ggtitle("Difference in dropout rates per cell") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom',axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank())
    
    library(foreach)
    Hicks_line_meanzeros<- foreach(i=1:length(sces),.combine=rbind)%do%{
        yy<-features$MeanRankZerosDiff[features$Dataset==names(sces)[i]]
        xx<-features$RankCounts[features$Dataset==names(sces)[i]]
        #plot(xx,yy,log='x')
        smDat <- sm(x = xx, y = yy,n.bins =25)
        qq<-cbind(smDat,rep(names(sces)[i],dim(smDat)[1]))
        #lines(smDat$x, smDat$y, lwd = 3, col = 4, lty=2)
        return(qq)
    }
    colnames(Hicks_line_meanzeros)[5]<-'Dataset'
    
    
    # yy<-features$MeanRankZerosDiff[features$Dataset==names(sces)[2]]
    # xx<-features$RankCounts[features$Dataset==names(sces)[2]]
    # plot(xx,yy)
    
    
    mean.zeros <- ggplot(features,aes_string(x = "RankCounts", y = "MeanRankZerosDiff", colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "black",alpha=hline_alpha,size=hline_size,linetype='dotted') +
        #geom_point(size = point.size, alpha = point.alpha,shape=46) +
        geom_smooth(size = linewidth)+
        #geom_line(data=Hicks_line_meanzeros,aes_string(x=Hicks_line_meanzeros$x,y=Hicks_line_meanzeros$y,colour = "Dataset"),size = linewidth.exp,alpha = 0.8)+
        scale_colour_manual(values = colours[-1]) +
        scale_fill_manual(values = colours[-1]) +
        xlab("Expression rank") +
        ylab("Difference in dropout rates") +
        ggtitle("Difference in mean-dropout rates relationship") +
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position ='bottom')
    
    return(list(mean.var=mean.var,z.gene=z.gene,z.cell=z.cell,mean.zeros=mean.zeros))
    #return(list(z.gene=z.gene,z.cell=z.cell))
}#end
