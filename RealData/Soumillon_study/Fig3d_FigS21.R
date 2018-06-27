#######DE fun#######
library(ROCR)
DE_smallfun<-function(LIST,TRUELABEL){
    
    auc_list<-list()
    for(i in 1:length(LIST)){
        
        auc_vec<-NULL
        
        for(j in 1:length(LIST[[i]])){
            
            if(length(dim(LIST[[i]][[j]])==2)){
                adjP<-apply(LIST[[i]][[j]],1,median)
                
                pred_MAST <- prediction(adjP, TRUELABEL)
                perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
                auc_temp<-performance( pred_MAST, measure='auc' )
                auc_temp<-auc_temp@y.values[[1]]
                auc_vec<-c(auc_vec,auc_temp)
                
                
                
            }  else{
                
                pred_MAST <- prediction(LIST[[i]][[j]]$adjpval, TRUELABEL)
                perf_MAST <- performance( pred_MAST, "tpr", "fpr" )
                
                auc_temp<-performance( pred_MAST, measure='auc' )
                auc_temp<-auc_temp@y.values[[1]]
                auc_vec<-c(auc_vec,auc_temp)
                
                
                #recall
                # recall_temp<-performance( pred_MAST, measure='acc', x.measure="cutoff" )
                # recall_temp@x.values[[1]][       which.max(recall_temp@"y.values"[[1]])]
                
                # qq<-LIST[[i]][[j]]$adjpval
                # qq[qq<0.05]=0
                # qq[qq>0]=1
                # xTab <-table(qq,TRUELABEL)
                
                
            }
        }
        auc_list[[i]]<-auc_vec
        
        
        
    }
    names(auc_list)<-names(LIST) 
    return(auc_list)
    
}

###begin analysis#########

load("E:/RNAseqProject/Soumillon_2014/D3T_smallsamples.RData")
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

boxplot(list(g1=aBay_out$BETA[[1]],g2=aBay_out$BETA[[2]]))




###bayNorm#######
load("E:/RNAseqProject/Soumillon_2014/smallgroup_bayNorm.RData")

DE_bay<-DE_smallfun(LIST=LIST_bayNorm,TRUELABEL=TRUE_LABEL_input)
DE_bay
unlist(lapply(DE_bay,mean))


######SAVER#########
#load("E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER_BETA.RData")



DE_SAVER<-DE_smallfun(LIST=LIST_SAVER,TRUELABEL=TRUE_LABEL_input)
DE_SAVER
unlist(lapply(DE_SAVER,mean))



####SCnorm#######
load("E:/RNAseqProject/Soumillon_2014/smallgroup_SCnorm.RData")

DE_SCnorm<-DE_smallfun(LIST=LIST_SCnorm,TRUELABEL=TRUE_LABEL_input)
DE_SCnorm
unlist(lapply(DE_SCnorm,mean))


####Scaling#######
load("E:/RNAseqProject/Soumillon_2014/smallgroup_Scaling.RData")
DE_Scaling<-DE_smallfun(LIST=LIST_Scaling,TRUELABEL=TRUE_LABEL_input)
DE_Scaling
unlist(lapply(DE_Scaling,mean))



#####scImpute#########
load("E:/RNAseqProject/Soumillon_2014/smallgroup_scImpute.RData")
DE_scImpute<-DE_smallfun(LIST=LIST_scImpute,TRUELABEL=TRUE_LABEL_input)
DE_scImpute

###MAGIC#########
load("E:/RNAseqProject/Soumillon_2014/smallgroup_MAGIC.RData")
DE_MAGIC<-DE_smallfun(LIST=LIST_MAGIC,TRUELABEL=TRUE_LABEL_input)
DE_MAGIC
unlist(lapply(DE_MAGIC,mean))

save(DE_bay,DE_SAVER,DE_SCnorm,DE_Scaling,DE_scImpute,DE_MAGIC,file="E:/RNAseqProject/Soumillon_2014/AUC_results.RData")

#######AUC_results v1 100 200 400############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_bayNorm.RData")
#load("E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER.RData")
#load("E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER_BETA.RData")
load("E:/RNAseqProject/Soumillon_2014/saver_default/smallgroup_SAVER_default.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_SCnorm.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_Scaling.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_scImpute.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_MAGIC.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_DCA.RData")



DE_bay<-DE_smallfun(LIST=LIST_bayNorm,TRUELABEL=TRUE_LABEL_input)
DE_SAVER<-DE_smallfun(LIST=LIST_SAVER,TRUELABEL=TRUE_LABEL_input)
DE_SCnorm<-DE_smallfun(LIST=LIST_SCnorm,TRUELABEL=TRUE_LABEL_input)
DE_Scaling<-DE_smallfun(LIST=LIST_Scaling,TRUELABEL=TRUE_LABEL_input)
DE_MAGIC<-DE_smallfun(LIST=LIST_MAGIC,TRUELABEL=TRUE_LABEL_input)
DE_scImpute<-DE_smallfun(LIST=LIST_scImpute,TRUELABEL=TRUE_LABEL_input)
DE_DCA<-DE_smallfun(LIST=LIST_DCA,TRUELABEL=TRUE_LABEL_input)

save(DE_bay,DE_SAVER,DE_SCnorm,DE_Scaling,DE_scImpute,DE_MAGIC,DE_DCA,file="E:/RNAseqProject/Soumillon_2014/AUC_results_default.RData")




###AUC_results 20 50 80#####
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_SCnorm_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_MAGIC_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_DCA_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_Scaling_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_scImpute_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/small_group_result_V2/smallgroup_bayNorm_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/saver_default/smallgroup_SAVER_default_v2.RData")

DE_bayNorm_v2<-DE_smallfun(LIST=LIST_bayNorm,TRUELABEL=TRUE_LABEL_input)
DE_SAVER_v2<-DE_smallfun(LIST=LIST_SAVER,TRUELABEL=TRUE_LABEL_input)
DE_SCnorm_v2<-DE_smallfun(LIST=LIST_SCnorm,TRUELABEL=TRUE_LABEL_input)
DE_scImpute_v2<-DE_smallfun(LIST=LIST_scImpute,TRUELABEL=TRUE_LABEL_input)
DE_Scaling_v2<-DE_smallfun(LIST=LIST_Scaling,TRUELABEL=TRUE_LABEL_input)
DE_MAGIC_v2<-DE_smallfun(LIST=LIST_MAGIC,TRUELABEL=TRUE_LABEL_input)
DE_DCA_v2<-DE_smallfun(LIST=LIST_DCA,TRUELABEL=TRUE_LABEL_input)

save(DE_bayNorm_v2,DE_SAVER_v2,DE_SCnorm_v2,DE_Scaling_v2,DE_scImpute_v2,DE_MAGIC_v2,DE_DCA_v2,file="E:/RNAseqProject/Soumillon_2014/AUC_results_v2.RData")



######begin plotting ##########
load("E:/RNAseqProject/Soumillon_2014/AUC_results_default.RData")
#load("E:/RNAseqProject/Soumillon_2014/AUC_results_1overBETA.RData")
#load("E:/RNAseqProject/Soumillon_2014/AUC_results_BETA.RData")
comparr<-c('100 vs 100','200 vs 200','400 vs 400','100 vs 200', '200 vs 100', '100 vs 400', '400 vs 100', '200 vs 400', '400 vs 200')


library(foreach)

method_vec<-c('bayNorm','SCnorm','scImpute','Scaling','SAVER','MAGIC','DCA')
AUC_list<-list(DE_bay,DE_SCnorm,DE_scImpute,DE_Scaling,DE_SAVER,DE_MAGIC,DE_DCA)
names(AUC_list)<-method_vec


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
col_vec<-cbbPalette[method_vec]

library(foreach)
boxdata<-foreach(i=1:length(AUC_list),.combine=rbind)%:%
    foreach(j=1:9,.combine=rbind)%do%{
        
        qqq<-cbind(AUC_list[[i]][[j]],rep(method_vec[i],length(AUC_list[[i]][[j]])),rep(comparr[j],length(AUC_list[[i]][[j]])))
        
        
    }


unlist(lapply(DE_bay,mean))
unlist(lapply(DE_MAGIC,mean))


boxdata<-as.data.frame(boxdata)
colnames(boxdata)<-c('AUC','Normalization methods','Group 1 vs Group 2')
boxdata[,1]<-as.numeric(as.character(boxdata[,1]))
boxdata[,3]<-factor(boxdata[,3],levels=unique(boxdata[,3]))
boxdata[,2]<-factor(boxdata[,2],levels=unique(boxdata[,2]))

cbbPalette_v3 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",'#CAB2D6','#6A3D9A','#FFFF99','#B15928')

library(ggplot2)
textsize<-10
BOX_PLOT<-ggplot(data=boxdata, aes(x=boxdata[,2], y=boxdata[,1], fill=boxdata[,3])) +
    geom_boxplot(size=0.05,outlier.size = 0.05)+
    #geom_line()+
    labs(fill = "# of cells: Group 1 vs Group 2",y='AUC',x='')+ggtitle("") +
    geom_hline(yintercept = 0.5,linetype=2) + 
    scale_colour_manual(values=col_vec )+
    #scale_colour_manual( values = cbbPalette_v3)+
    #guides (fill = guide_legend (override.aes=list(fill=NA)))+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top',legend.background=element_blank(),legend.key =element_blank())

BOX_PLOT


BOX_PLOT
ggsave(filename="E:/RNAseqProject/Illustrator_bayNorm/FIGURE_2/Soumillon_smallgroups_default.pdf",units='in',plot=BOX_PLOT,width=7,height=4,device='pdf')



#######Fig 3 (d)##########
load("E:/RNAseqProject/Soumillon_2014/AUC_results_default.RData")
load("E:/RNAseqProject/Soumillon_2014/AUC_results_v2.RData")

comparr<-c('20 vs 20','50 vs 50','80 vs 80','100 vs 100','200 vs 200','400 vs 400')


library(foreach)

method_vec<-c('bayNorm','SCnorm','scImpute','Scaling','SAVER','MAGIC','DCA')


AUC_list<-list(c(DE_bayNorm_v2[seq(1:3)],DE_bay[seq(1:3)]),
               c(DE_SCnorm_v2[seq(1:3)],DE_SCnorm[seq(1:3)]),
               c(DE_scImpute_v2[seq(1:3)],DE_scImpute[seq(1:3)]),
               c(DE_Scaling_v2[seq(1:3)],DE_Scaling[seq(1:3)]),
               c( DE_SAVER_v2[seq(1:3)] ,DE_SAVER[seq(1:3)]),
               c(DE_MAGIC_v2[seq(1:3)],DE_MAGIC[seq(1:3)]),
               c(DE_DCA_v2[seq(1:3)],DE_DCA[seq(1:3)]))
names(AUC_list)<-method_vec


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
col_vec<-cbbPalette[method_vec]



library(foreach)
boxdata<-foreach(i=1:length(AUC_list),.combine=rbind)%:%
    foreach(j=1:length(comparr),.combine=rbind)%do%{
        
        qqq<-cbind(AUC_list[[i]][[j]],rep(method_vec[i],length(AUC_list[[i]][[j]])),rep(comparr[j],length(AUC_list[[i]][[j]])))
        
        
    }


unlist(lapply(DE_bay,mean))
unlist(lapply(DE_MAGIC,mean))


boxdata<-as.data.frame(boxdata)
colnames(boxdata)<-c('AUC','Normalization methods','Group 1 vs Group 2')
boxdata[,1]<-as.numeric(as.character(boxdata[,1]))
boxdata[,3]<-factor(boxdata[,3],levels=unique(boxdata[,3]))
boxdata[,2]<-factor(boxdata[,2],levels=unique(boxdata[,2]))

cbbPalette_v3 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",'#CAB2D6','#6A3D9A','#FFFF99','#B15928')

library(ggplot2)
textsize<-10
BOX_PLOT<-ggplot(data=boxdata, aes(x=boxdata[,2], y=boxdata[,1], fill=boxdata[,3])) +
    geom_boxplot(size=0.05,outlier.size = 0.05)+
    #facet_wrap(~as.factor(boxdata[,3]), nrow=1)+
    #geom_line()+
    labs(fill = "# of cells: Group 1 vs Group 2",y='AUC',x='')+ggtitle("") +
    geom_hline(yintercept = 0.5,linetype=2) +
    #scale_fill_manual(values=col_vec )+
    scale_fill_manual( values = cbbPalette_v3)+
    scale_x_discrete(limits=c('bayNorm','DCA','MAGIC','SAVER','SCnorm','scImpute','Scaling'))+
    #guides (fill = guide_legend (override.aes=list(fill=NA)))+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top',legend.background=element_blank(),legend.key =element_blank())

BOX_PLOT


ggsave(filename="E:/RNAseqProject/Illustrator_bayNorm/FIGURE_2/Soumillon_smallgroups_default_main.pdf",units='in',plot=BOX_PLOT,width=7,height=4,device='pdf')





#########SUP Fig S21##############
comparr<-c('20 vs 50','50 vs 20','20 vs 80', '80 vs 20', '50 vs 80','80 vs 50','100 vs 200', '200 vs 100', '100 vs 400', '400 vs 100', '200 vs 400', '400 vs 200')


library(foreach)

method_vec<-c('bayNorm','SCnorm','scImpute','Scaling','SAVER','MAGIC','DCA')


AUC_list<-list(c(DE_bayNorm_v2[-seq(1:3)],DE_bay[-seq(1:3)]),
               c(DE_SCnorm_v2[-seq(1:3)],DE_SCnorm[-seq(1:3)]),
               c(DE_scImpute_v2[-seq(1:3)],DE_scImpute[-seq(1:3)]),
               c(DE_Scaling_v2[-seq(1:3)],DE_Scaling[-seq(1:3)]),
               c( DE_SAVER_v2[-seq(1:3)] ,DE_SAVER[-seq(1:3)]),
               c(DE_MAGIC_v2[-seq(1:3)],DE_MAGIC[-seq(1:3)]),
               c(DE_DCA_v2[-seq(1:3)],DE_DCA[-seq(1:3)]))
names(AUC_list)<-method_vec


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
col_vec<-cbbPalette[method_vec]

library(foreach)
boxdata<-foreach(i=1:length(AUC_list),.combine=rbind)%:%
    foreach(j=1:length(comparr),.combine=rbind)%do%{
        
        qqq<-cbind(AUC_list[[i]][[j]],rep(method_vec[i],length(AUC_list[[i]][[j]])),rep(comparr[j],length(AUC_list[[i]][[j]])))
        
        
    }


unlist(lapply(DE_bay,mean))
unlist(lapply(DE_MAGIC,mean))


boxdata<-as.data.frame(boxdata)
colnames(boxdata)<-c('AUC','Normalization methods','Group 1 vs Group 2')
boxdata[,1]<-as.numeric(as.character(boxdata[,1]))
boxdata[,3]<-factor(boxdata[,3],levels=unique(boxdata[,3]))
boxdata[,2]<-factor(boxdata[,2],levels=unique(boxdata[,2]))

cbbPalette_v3 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",'#CAB2D6','#6A3D9A','#FFFF99','#B15928')

library(ggplot2)
textsize<-10
BOX_PLOT<-ggplot(data=boxdata, aes(x=boxdata[,2], y=boxdata[,1], fill=boxdata[,3])) +
    geom_boxplot(size=0.05,outlier.size = 0.05)+
    #geom_line()+
    labs(fill = "# of cells: Group 1 vs Group 2",y='AUC',x='')+ggtitle("") +
    geom_hline(yintercept = 0.5,linetype=2) + 
    #scale_fill_manual(values=col_vec )+
    scale_fill_manual( values = cbbPalette_v3)+
    scale_x_discrete(limits=c('bayNorm','DCA','MAGIC','SAVER','SCnorm','scImpute','Scaling'))+
    #guides (fill = guide_legend (override.aes=list(fill=NA)))+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top',legend.background=element_blank(),legend.key =element_blank())

BOX_PLOT
ggsave(filename="E:/RNAseqProject/Illustrator_bayNorm/FIGURE_2/Soumillon_smallgroups_default_sup.pdf",units='in',plot=BOX_PLOT,width=8,height=6,device='pdf')
