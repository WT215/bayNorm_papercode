###CV functions######

CV_fun <- function(x){  
    qq<-sd(x) / mean(x) 
    return(qq)
}
bootstrap_fun <- function(inputdata,smFISH_genename, type = '') {
    boot_mat <- matrix(nrow = 1000, ncol = length(smFISH_genename))
    colnames(boot_mat) <- smFISH_genename
    
    library(foreach)
    library(doSNOW)
    library(parallel)
    
    cluster = makeCluster(5, type = "SOCK")
    registerDoSNOW(cluster)
    getDoParWorkers()
    
    
    if (type == 'smFISH') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c) %do% {
                qq_temp <- CV_fun(sample(smFISH_list[[which(names(smFISH_list)==smFISH_genename[i])]], replace = T))
                return(qq_temp)
            }
        }
    } else if (type == 'array3D') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c,.export='CV_fun') %dopar% {
                qq_temp <- CV_fun(sample(as.numeric(inputdata[smFISH_genename[i], , ]), replace = T))
                return(qq_temp)
            }
        }
        
    } else if (type == 'mat2D') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c,.export='CV_fun') %dopar% {
                qq_temp <- CV_fun(sample(inputdata[smFISH_genename[i], ], replace = T))
                return(qq_temp)
            }
        }
        
    }
    stopCluster(cluster)
    
    
    
    return(boot_mat)
    
}


####density functions#######
qwerfun_LMNA<-function(DenDat,linesize=1,CAPTION=''){
    
    pp<-unique(DenDat$`Normalization method`)
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbbPalette )<-c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','Raw')
    cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% pp)]
    
    
    Gg<-unique(DenDat$Genename)
    typp<-unique(DenDat$Condition)
    
    filenn<-paste('Gene: ',Gg,sep='')
    
    density_plot<-ggplot() +
        stat_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,geom="line",position="identity")+
        scale_colour_manual( values = cbbPalette2)+
        #add xlim for other than Grun data
        xlim(0, 500)+
        ylim(0,0.008)+
        labs(x = "Gene expression",y='density',color='Normalization methods',caption=CAPTION)+
        ggtitle(filenn)+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position = c(0.8,0.9))
    return(density_plot)
}


qwerfun_ind_LMNA<-function(DenDat,filenn='',linesize=1,CAPTION=''){
    
    pp<-unique(DenDat$`Normalization method`)
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbbPalette )<-c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','Raw')
    cbbPalette2<-cbbPalette[match(pp,names(cbbPalette))]
    
    Gg<-unique(DenDat$Genename)
    typp<-unique(DenDat$Condition)
    
    filenn<-paste('Gene: ',Gg,sep='')
    
    density_plot<-ggplot() +
        stat_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,geom="line",position="identity",show.legend = F)+
        xlim(0, 500)+
        ylim(0,0.018)+
        scale_colour_manual( values = cbbPalette2)+
        labs(x = "Normalized count",y='density',color='Normalization methods')+
        ggtitle(filenn,subtitle=CAPTION)+
        guides (colour = guide_legend (override.aes=list(colour=NA)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1.5,"line"),legend.position = c(0.85,0.85),legend.background=element_blank(),legend.key =element_blank())
    
    g1 <- ggplotGrob(density_plot)
    # Check out the grobs
    library(grid)
    #grid.ls(grid.force(g1))
    # Get names of 'label' grobs.
    names.grobs <- grid.ls(grid.force(g1), print=F)$name 
    labels <- names.grobs[which(grepl("^label", names.grobs))]
    
    for(i in seq_along(labels)) {
        g1 <- editGrob(grid.force(g1), gPath(labels[i]), grep = TRUE,  gp = gpar(col =cbbPalette2[i]))
    }
    
    # Draw it
    #grid.newpage()
    #density_plot_v2<- grid.draw(g1, recording=F)
    density_plot_v2<-as_ggplot(g1)
    
    
    return(density_plot_v2)
}


######begin analysis###########

load("Torre_many_normalizations.RData")

#scale normalized data by mean BETA used in bayNorm, for a fair comparison####
MAGIC_out<-as.matrix(MAGIC_Torre)
MAGIC_out<-MAGIC_out/mean(bay_out$BETA)
scImpute_out<-scImpute_out/mean(bay_out$BETA)
SCnorm_dat<-SCnorm_out@metadata$NormalizedData/mean(bay_out$BETA)
saver_array<-saver_array/mean(bay_out$BETA)

DCA_Torre<-DCA_Torre/mean(bay_out$BETA)


bayinput<-bay_out$Bay_array

smFISH_genename<-intersect(rownames(bayinput),colnames(Torre_FISH_sub_norm))
length(smFISH_genename)

Raw_normal<-Torre_drop_sub
RB_norm <- t(t(Torre_drop_sub) / bay_out$BETA)



#begin bootstrapping (bootstrapping results are not shown in the paper)####
library(foreach)
set.seed(12300)
boot_smFISH <- bootstrap_fun(inputdata=NULL, smFISH_genename,type = 'smFISH')
boot_bayNorm <-bootstrap_fun(inputdata = bay_out$Bay_array[smFISH_genename,,], smFISH_genename,type = 'array3D')
boot_saver <-bootstrap_fun(inputdata =saver_array[smFISH_genename,,],smFISH_genename, type = 'array3D')
boot_RB <- bootstrap_fun(inputdata = RB_norm[smFISH_genename,],smFISH_genename, type = 'mat2D')
#boot_raw <- bootstrap_fun(inputdata = Raw_normal,smFISH_genename, type = 'mat2D')
boot_scnorm <- bootstrap_fun(inputdata =SCnorm_dat[smFISH_genename,],smFISH_genename, type = 'mat2D')
boot_scImpute <- bootstrap_fun(inputdata =scImpute_out[smFISH_genename,],smFISH_genename, type = 'mat2D')
boot_magic <- bootstrap_fun(inputdata =MAGIC_out[smFISH_genename,],smFISH_genename, type = 'mat2D')
boot_dca <- bootstrap_fun(inputdata =DCA_Torre[smFISH_genename,],smFISH_genename, type = 'mat2D')

save(boot_smFISH,boot_bayNorm,boot_saver,boot_RB,boot_scnorm,boot_scImpute,boot_magic, boot_dca ,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/boot_results_tr_default_divmeanbeta.RData")




load("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/boot_results_tr_default_divmeanbeta.RData")

#####CV MEAN ERROR#####
CV_smFISH <-
    unlist(lapply(smFISH_list, function(x) {
        CV_fun(x)
    }))[smFISH_genename]
ERROR_smFISH <- 1.96 * sqrt(apply(boot_smFISH, 2, var))[smFISH_genename]
MEAN_smFISH <- colMeans(boot_smFISH)[smFISH_genename]

Mu_Test <- apply(bay_out$Bay_array[smFISH_genename, , ], 1, mean)
Sd_Test <-
    apply(apply(bay_out$Bay_array[smFISH_genename, , ], c(1, 3), sd), 1, mean)
CV_Test <- Sd_Test / Mu_Test
ERROR_Test <- 1.96 * sqrt(apply(boot_bayNorm, 2, var))[smFISH_genename]
MEAN_Test <- colMeans(boot_bayNorm)[smFISH_genename]

Mu_saver <- apply(saver_array[smFISH_genename, , ], 1, mean)
Sd_saver <-
    apply(apply(saver_array[smFISH_genename, , ], c(1, 3), sd), 1, mean)
CV_saver <- Sd_saver / Mu_saver
ERROR_saver <- 1.96 * sqrt(apply(boot_saver, 2, var))[smFISH_genename]
MEAN_saver <- colMeans(boot_saver)[smFISH_genename]



CV_RB <- apply(RB_norm[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_RB <- 1.96 * sqrt(apply(boot_RB, 2, var))[smFISH_genename]
MEAN_RB <- colMeans(boot_RB)[smFISH_genename]

CV_scnorm <- apply(SCnorm_dat[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_scnorm <- 1.96 * sqrt(apply(boot_scnorm, 2, var))[smFISH_genename]
MEAN_scnorm <- colMeans(boot_scnorm)[smFISH_genename]

CV_scnorm <- apply(SCnorm_out@metadata$NormalizedData[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_scnorm <- 1.96 * sqrt(apply(boot_scnorm, 2, var))[smFISH_genename]
MEAN_scnorm <- colMeans(boot_scnorm)[smFISH_genename]


CV_scimpute <- apply(scImpute_out[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_scimpute <- 1.96 * sqrt(apply(boot_scImpute, 2, var))[smFISH_genename]
MEAN_scimpute <- colMeans(boot_scImpute)[smFISH_genename]

CV_magic <- apply(MAGIC_out[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_magic <- 1.96 * sqrt(apply(boot_magic, 2, var))[smFISH_genename]
MEAN_magic<- colMeans(boot_magic)[smFISH_genename]

CV_dca <- apply(DCA_Torre[smFISH_genename, ], 1, function(x) {
    CV_fun(x)
})
ERROR_dca <- 1.96 * sqrt(apply(boot_dca, 2, var))[smFISH_genename]
MEAN_dca<- colMeans(boot_dca)[smFISH_genename]


CV_vec <- c(CV_smFISH, CV_Test,CV_scnorm,CV_RB,CV_saver,CV_scimpute,CV_magic,CV_dca)
ERROR_vec <- c(ERROR_smFISH, ERROR_Test,ERROR_scnorm,ERROR_RB,ERROR_saver,ERROR_scimpute,ERROR_magic,ERROR_dca)
MEAN_vec <- c(MEAN_smFISH, MEAN_Test,MEAN_scnorm,MEAN_RB,MEAN_saver,MEAN_scimpute,MEAN_magic,MEAN_dca)


lg<-length(smFISH_genename)
Methods <- c(rep('smFISH',lg), rep('bayNorm',lg),rep('SCnorm',lg),rep('Scaling',lg), rep('SAVER',lg),rep('scImpute',lg),rep('MAGIC',lg),rep('DCA',lg))


BAR_DAT <-data.frame(
        CV = CV_vec,
        Method = Methods,
        Genename = rep(smFISH_genename, length(unique(Methods))),
        ERROR = ERROR_vec,
        MEAN = MEAN_vec
    )
CAPTION <- 'Single cell data based on Torre case study'


BAR_DAT_input <- BAR_DAT
BAR_DAT_input[, 1] <- as.numeric(as.character(BAR_DAT_input[, 1]))
BAR_DAT_input[, 3] <-
    factor(BAR_DAT_input[, 3], levels = unique(BAR_DAT_input[, 3]))
BAR_DAT_input[, 2] <-
    factor(BAR_DAT_input[, 2], levels = unique(BAR_DAT_input[, 2]))



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
textsize<-12
BAR_out<-ggplot(data=BAR_DAT_input, aes(x=BAR_DAT_input[,3], y=2*BAR_DAT_input$CV-BAR_DAT_input$MEAN, fill=BAR_DAT_input[,2])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_errorbar(aes(ymin=2*BAR_DAT_input$CV-BAR_DAT_input$MEAN-BAR_DAT_input$ERROR, ymax=2*BAR_DAT_input$CV-BAR_DAT_input$MEAN+BAR_DAT_input$ERROR),size=.5,width=.4,position=position_dodge(.9))+
    labs(x = "Gene",y='CV',fill='Methods',caption=CAPTION)+
    scale_fill_manual(values=cbbPalette)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
BAR_out




#######make density plots############
library(foreach)
DAT_LIST<-list(smFISH=smFISH_list,bayNorm=bay_out$Bay_array[smFISH_genename,,],SAVER=saver_array[smFISH_genename,,],Scaling=RB_norm[smFISH_genename,],Raw=Torre_drop_sub[smFISH_genename,],MAGIC=MAGIC_out[smFISH_genename,],scImpute=scImpute_out[smFISH_genename,],SCnorm=SCnorm_out@metadata$NormalizedData[smFISH_genename,],DCA=DCA_Torre[smFISH_genename,])

CNF<-function(x,numb=500){
    qq<-x
    return(qq)
}

DenDat_list<-list()
for(i in 1:length(smFISH_genename)){
    print(i)
    Gg<-smFISH_genename[i]
    
    input_smfish<-smFISH_list[[which(names(smFISH_list)==Gg)]]
    input_bay<-as.numeric(DAT_LIST$bayNorm[Gg,,])
    input_SAVER<-as.numeric(DAT_LIST$SAVER[Gg,,])
    input_Scaling=RB_norm[Gg,]
    input_raw<-DAT_LIST$Raw[Gg,]
    input_scImpute<-DAT_LIST$scImpute[Gg,]
    input_MAGIC<-DAT_LIST$MAGIC[Gg,]
    input_SCnorm<-DAT_LIST$SCnorm[Gg,]
    input_DCA<-DAT_LIST$DCA[Gg,]
    
    inputlist<-list(smFISH=input_smfish,bayNorm=input_bay,SAVER=input_SAVER,Scaling=input_Scaling,Raw=input_raw)
    #inputlist<-list(smFISH=input_smfish,bayNorm=input_bay,Scaling=input_Scaling,Raw=input_raw)
    #inputlist<-list(smFISH=input_smfish,bayNorm=input_bay,SAVER=input_SAVER,Scaling=input_Scaling,Raw=input_raw,scImpute=input_scImpute,MAGIC=input_MAGIC,SCnorm=input_SCnorm)
    methodd<-names(inputlist)
    
    library(foreach)
    DenDat<-foreach(kkk=1:length(inputlist),.combine=rbind)%do%{
        aa<-inputlist[[kkk]]
        bb<-cbind(aa,rep(methodd[kkk],length(aa)),rep(kkk,length(aa)),rep(Gg,length(aa)))
        return(bb)
    }
    colnames(DenDat)<-c('Normalized count','Normalization method','Colour','Genename')
    DenDat<-as.data.frame(DenDat)
    DenDat$`Normalized count`<-as.numeric(as.character(DenDat$`Normalized count`))
    DenDat$`Normalization method`<-factor(DenDat$`Normalization method`,levels=unique(DenDat$`Normalization method`))
    
    
    DenDat_list[[i]]<-DenDat
}


source("E:/RNAseqProject/Density_fun.r")
library(gridExtra)
library(ggpubr)
bw<-10
textsize<-10
#cbPalette <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gg_list<-lapply(DenDat_list,qwerfun,CAPTION= '',linesize=0.5)
names(gg_list)<-smFISH_genename
names( DenDat_list)<-smFISH_genename

######multiple density plots#########

qq<-do.call('ggarrange',c(gg_list,ncol=3,nrow=4,common.legend = TRUE, legend="bottom"))
qq
ggsave(filename ="E:/RNAseqProject/Illustrator_bayNorm/SUP/smFISH/Torre_den_default_divmeanbeta.pdf",plot = qq,width=8,height=9)


####single LMNA#########

textsize<-12
#linesize=1 for individual LMNA plot
gg_list<-lapply(DenDat_list,qwerfun,CAPTION= '',linesize=1)
names(gg_list)<-smFISH_genename
names( DenDat_list)<-smFISH_genename

LMNA_plot<-qwerfun_ind_LMNA(DenDat_list$LMNA,CAPTION='Data from Torre et al study')
LMNA_plot


###The following code needs to be run for making Figures 2(c)-(h)


#########mean#######
MeanExp_smFISH<-cbind(unlist(lapply(smFISH_list,mean))[smFISH_genename],smFISH_genename,rep('smFISH',length(smFISH_genename)))
MeanExp_bayNorm<-cbind(rowMeans(bay_out$Bay_array[smFISH_genename,,])[smFISH_genename],smFISH_genename,rep('bayNorm',length(smFISH_genename)))
MeanExp_saver<-cbind(rowMeans(saver_array[smFISH_genename,,])[smFISH_genename],smFISH_genename,rep('SAVER',length(smFISH_genename)))
MeanExp_SCnorm<-cbind(rowMeans(SCnorm_dat[smFISH_genename,]),smFISH_genename,rep('SCnorm',length(smFISH_genename)))
MeanExp_Scaling<-cbind(rowMeans(RB_norm[smFISH_genename,]),smFISH_genename,rep('Scaling',length(smFISH_genename)))
MeanExp_scImpute<-cbind(rowMeans(scImpute_out[smFISH_genename,]),smFISH_genename,rep('scImpute',length(smFISH_genename)))
MeanExp_MAGIC<-cbind(rowMeans(MAGIC_out[smFISH_genename,]),smFISH_genename,rep('MAGIC',length(smFISH_genename)))
MeanExp_DCA<-cbind(rowMeans(DCA_Torre[smFISH_genename,]),smFISH_genename,rep('DCA',length(smFISH_genename)))



MeanExp_data_Torre<-rbind(MeanExp_smFISH,MeanExp_bayNorm,MeanExp_saver,MeanExp_SCnorm,MeanExp_Scaling,MeanExp_scImpute,MeanExp_MAGIC,MeanExp_DCA)
colnames(MeanExp_data_Torre)<-c('MeanExp','Gene','NormMethods')


save(MeanExp_data_Torre,BAR_DAT,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/Torre_CV_default_divmeanbeta.RData")

