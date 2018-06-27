bootstrap_fun <- function(inputdata,smFISH_genename, type = '') {
    boot_mat <- matrix(nrow = 1000, ncol = length(smFISH_genename))
    colnames(boot_mat) <- smFISH_genename
    library(ineq)
    
    library(foreach)
    library(doSNOW)
    library(parallel)
    
    cluster = makeCluster(5, type = "SOCK")
    registerDoSNOW(cluster)
    getDoParWorkers()
    
    
    if (type == 'smFISH') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c) %do% {
                qq_temp <- Gini(sample(smFISH_list[[which(names(smFISH_list)==smFISH_genename[i])]], replace = T))
                return(qq_temp)
            }
        }
    } else if (type == 'array3D') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c,.packages = 'ineq') %dopar% {
                qq_temp <- Gini(sample(as.numeric(inputdata[smFISH_genename[i], , ]), replace = T))
                return(qq_temp)
            }
        }
        
    } else if (type == 'mat2D') {
        for (i in 1:dim(boot_mat)[2]) {
            boot_mat[, i] <- foreach(qq = 1:1000, .combine = c,.packages = 'ineq') %dopar% {
                qq_temp <- Gini(sample(inputdata[smFISH_genename[i], ], replace = T))
                return(qq_temp)
            }
        }
        
    }
    stopCluster(cluster)
    
    
    
    return(boot_mat)
    
}



load("Torre_many_normalizations.RData")


MAGIC_out<-as.matrix(MAGIC_Torre)
bayinput<-bay_out$Bay_array
smFISH_genename<-intersect(rownames(bayinput),colnames(Torre_FISH_sub_norm))
length(smFISH_genename)

Raw_normal<-Torre_drop_sub
#Scaling method
RB_norm <- t(t(Torre_drop_sub) / bay_out$BETA)

MAGIC_out<-as.matrix(MAGIC_Torre)
MAGIC_out<-MAGIC_out/mean(bay_out$BETA)
scImpute_out<-scImpute_out/mean(bay_out$BETA)
SCnorm_dat<-SCnorm_out@metadata$NormalizedData/mean(bay_out$BETA)
saver_array<-saver_array/mean(bay_out$BETA)
#DCA was implemented in python, so code for DCA was provided in other files
DCA_Torre<-DCA_Torre/mean(bay_out$BETA)


###bootstrapping: estimate Gini####
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

save(boot_smFISH,boot_bayNorm,boot_saver,boot_RB,boot_scnorm,boot_scImpute,boot_magic,boot_dca,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/boot_results_GINI_default.RData")

load("E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/boot_results_GINI_default.RData")



#####Gini MEAN ERROR#####
Gini_smFISH <-
    unlist(lapply(smFISH_list, function(x) {
        Gini(x)
    }))[smFISH_genename]
ERROR_smFISH <- 1.96 * sqrt(apply(boot_smFISH, 2, var))[smFISH_genename]
MEAN_smFISH <- colMeans(boot_smFISH)[smFISH_genename]

Gini_bay<-apply(bay_out$Bay_array[smFISH_genename,,],1,Gini)[smFISH_genename]
ERROR_Test <- 1.96 * sqrt(apply(boot_bayNorm, 2, var))[smFISH_genename]
MEAN_Test <- colMeans(boot_bayNorm)[smFISH_genename]

Gini_saver<-apply(saver_array[smFISH_genename,,],1,Gini)[smFISH_genename]
ERROR_saver <- 1.96 * sqrt(apply(boot_saver, 2, var))[smFISH_genename]
MEAN_saver <- colMeans(boot_saver)[smFISH_genename]



Gini_RB <- apply(RB_norm[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_RB <- 1.96 * sqrt(apply(boot_RB, 2, var))[smFISH_genename]
MEAN_RB <- colMeans(boot_RB)[smFISH_genename]

# Gini_raw <- apply(Raw_normal[smFISH_genename, ], 1, function(x) {
#   Gini(x)
# })
# ERROR_raw <- 1.96 * sqrt(apply(boot_raw, 2, var))[smFISH_genename]
# MEAN_raw <- colMeans(boot_raw)[smFISH_genename]

Gini_scnorm <- apply(SCnorm_out@metadata$NormalizedData[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_scnorm <- 1.96 * sqrt(apply(boot_scnorm, 2, var))[smFISH_genename]
MEAN_scnorm <- colMeans(boot_scnorm)[smFISH_genename]

Gini_scnorm <- apply(SCnorm_out@metadata$NormalizedData[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_scnorm <- 1.96 * sqrt(apply(boot_scnorm, 2, var))[smFISH_genename]
MEAN_scnorm <- colMeans(boot_scnorm)[smFISH_genename]


Gini_scimpute <- apply(scImpute_out[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_scimpute <- 1.96 * sqrt(apply(boot_scImpute, 2, var))[smFISH_genename]
MEAN_scimpute <- colMeans(boot_scImpute)[smFISH_genename]

Gini_magic <- apply(MAGIC_out[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_magic <- 1.96 * sqrt(apply(boot_magic, 2, var))[smFISH_genename]
MEAN_magic<- colMeans(boot_magic)[smFISH_genename]


Gini_dca<- apply(DCA_Torre[smFISH_genename, ], 1, function(x) {
    Gini(x)
})
ERROR_dca <- 1.96 * sqrt(apply(boot_dca, 2, var))[smFISH_genename]
MEAN_dca<- colMeans(boot_dca)[smFISH_genename]


Gini_vec <- c(Gini_smFISH, Gini_bay,Gini_scnorm,Gini_RB,Gini_saver,Gini_scimpute,Gini_magic,Gini_dca)
ERROR_vec <- c(ERROR_smFISH, ERROR_Test,ERROR_scnorm,ERROR_RB,ERROR_saver,ERROR_scimpute,ERROR_magic,ERROR_dca)
MEAN_vec <- c(MEAN_smFISH, MEAN_Test,MEAN_scnorm,MEAN_RB,MEAN_saver,MEAN_scimpute,MEAN_magic,MEAN_dca)

# CV_vec <- c(CV_smFISH, CV_Test, CV_saver,CV_RB, CV_raw)
# ERROR_vec <- c(ERROR_smFISH, ERROR_Test, ERROR_saver,ERROR_RB, ERROR_raw)
# MEAN_vec <- c(MEAN_smFISH, MEAN_Test, MEAN_saver,MEAN_RB, MEAN_raw)

lg<-length(smFISH_genename)
Methods <- c(rep('smFISH',lg), rep('bayNorm',lg),rep('SCnorm',lg),rep('Scaling',lg), rep('SAVER',lg),rep('scImpute',lg),rep('MAGIC',lg),rep('DCA',lg))


BAR_DAT <-
    data.frame(
        Gini = Gini_vec,
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
BAR_out<-ggplot(data=BAR_DAT_input, aes(x=BAR_DAT_input[,3], y=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN, fill=BAR_DAT_input[,2])) +
    #BAR_out<-ggplot(data=BAR_DAT_input, aes(x=BAR_DAT_input[,3], y=BAR_DAT_input$Gini, fill=BAR_DAT_input[,2])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    geom_errorbar(aes(ymin=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN-BAR_DAT_input$ERROR, ymax=2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN+BAR_DAT_input$ERROR),size=.5,width=.4,position=position_dodge(.9))+
    #geom_text(aes(label=round(2*BAR_DAT_input$Gini-BAR_DAT_input$MEAN,2)), vjust=1.6, color="black", position = position_dodge(0.9), size=2.5)+
    labs(x = "Gene",y='Gini',fill='Methods',caption=CAPTION)+
    #ggtitle('') +
    #scale_fill_brewer(palette="Paired")+
    scale_fill_manual(values=cbbPalette)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
#dev.off()

BAR_out


BAR_out_Torre_GINI<-BAR_DAT
#Save for making Figures 2 (g)-(h)
save( BAR_out_Torre_GINI,file="E:/RNAseqProject/NEWPROJECT_PAPERS/Torre_2017/BAY_8640_V2/New/Torre_GINI_default_divmeanbeta.RData")