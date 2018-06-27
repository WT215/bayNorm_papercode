source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

library(foreach)
M_TEST_FUN<-function(DATA_LIST,grs,methodnames=NULL){


    
    
  M_list<-foreach(i=1:length(DATA_LIST))%do%{

    print(i)
     if(length(dim(DATA_LIST[[i]]))==2){
         mmm<-SCnorm_runMAST(Data=DATA_LIST[[i]], NumCells=as.numeric(grs))
         cat("\014") 
         mmm_p<-mmm$adjpval
         return(mmm_p)
     } else if(length(dim(DATA_LIST[[i]]))==3){
         
         library(foreach)
         library(doSNOW)
         library(parallel)
         
         cluster = makeCluster(5, type = "SOCK")
         registerDoSNOW(cluster)
         getDoParWorkers()
         
         iterations <- dim(DATA_LIST[[i]])[3]
         pb <- txtProgressBar(max = iterations, style = 3)
         progress <- function(n) setTxtProgressBar(pb, n)
         opts <- list(progress = progress)
         
         wer<-foreach(sampleind=1:dim(DATA_LIST[[i]])[3],.combine=cbind,.options.snow = opts, .export=c('SCnorm_runMAST', 'SCnorm_runMAST2'), .packages=c('MAST','reshape2'))%dopar%{
             print(paste('sampleind',sampleind))
             mmm<-SCnorm_runMAST(Data=DATA_LIST[[i]][,,sampleind], NumCells=as.numeric(grs))
             cat("\014") 
             mmm_p<-mmm$adjpval
             return(mmm_p)
         }
         
         close(pb)
         stopCluster(cluster)
         
         return(apply(wer,1,median))
         
     }
      

    
    
  }

  
  
  names(M_list)<-methodnames
  return(M_list)

}


T_TEST_FUN<-function(DATA_LIST,CONDITION,methodnames=NULL){
  qqq<-unique(CONDITION)

  T_list<-foreach(i=1:length(DATA_LIST))%do%{

    print(i)
    # if(i==1){
    # 
    #   tttt<-tStatAnalysis_2groups_2(cells=DATA_LIST[[i]][,CONDITION==qqq[1],],ctrls=DATA_LIST[[i]][,CONDITION==qqq[2],], list_mode=NULL, verbose=T, plotFDR = FALSE)
    # 
    # }else{
      DATA_LIST[[i]]<-as.matrix(DATA_LIST[[i]])
      tttt<-tStatAnalysis_2groups_2(cells=DATA_LIST[[i]][,CONDITION==qqq[1]],ctrls=DATA_LIST[[i]][,CONDITION==qqq[2]], list_mode=NULL, verbose=T, plotFDR = FALSE)

   # }

    metric=5
    tttt_p<-do.call(cbind,lapply(tttt,function(x){return(x[[metric]])}))
    tttt_p2<-apply(tttt_p,1,median)
    return(tttt_p2)
  }

  names(T_list)<-methodnames
  return(T_list)

}



BATCH_BAR_FUN<-function(M_n1_list,methodnames=NULL){

method_vec<-methodnames
col_vec<-seq(1,length(method_vec))
library(foreach)
thrrr_vec<-c(0.01,0.05,0.1)
BAR_DAT<-foreach(i = seq(1,length(thrrr_vec)),.combine=rbind)%do%{
  thrrr<-thrrr_vec[i]
  ng_vec<-unlist(lapply(M_n1_list,function(x){length(which(x<thrrr))}))
  ng_vec<-cbind(ng_vec,rep(thrrr,length(M_n1_list)))
  ng_vec<-cbind(ng_vec,method_vec)
}



colnames(BAR_DAT)<-c('Number of detected DE genes','Adjusted P-values threshold','Normalization method')

BAR_DAT<-as.data.frame(BAR_DAT)

BAR_DAT[,1]<-as.numeric(as.character(BAR_DAT[,1]))
BAR_DAT[,2]<-factor(BAR_DAT[,2],levels=unique(BAR_DAT[,2]))
BAR_DAT[,3]<-factor(BAR_DAT[,3],levels=unique(BAR_DAT[,3]))

num_gene<-length(M_n1_list[[1]])

textsize<-40
N1_BATCH_FIG<-ggplot(data=BAR_DAT, aes(x=BAR_DAT[,2], y=BAR_DAT[,1]/num_gene, fill=BAR_DAT[,3])) +
  geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
  geom_text(aes(label=round(BAR_DAT[,1]/num_gene,4)), vjust=1.6, color="black", position = position_dodge(0.9), size=10)+
  labs(x = 'Adjusted P-values threshold',y='False positive rates',fill='Normalization methods')+ggtitle("(a)") +
  scale_fill_brewer(palette="Paired")+
  theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size =textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) )
return(N1_BATCH_FIG)
}




BATCH_BAR_FUN2<-function(M_n1_list,methodnames=NULL,textsize=6,geom_text_size=1.2,legend_key_size=0.5){
    
    method_vec<-methodnames
    col_vec<-seq(1,length(method_vec))
    library(foreach)
    thrrr_vec<-c(0.01,0.05,0.1)
    BAR_DAT<-foreach(i = seq(1,length(thrrr_vec)),.combine=rbind)%do%{
        thrrr<-thrrr_vec[i]
        ng_vec<-unlist(lapply(M_n1_list,function(x){length(which(x<thrrr))}))
        ng_vec<-cbind(ng_vec,rep(thrrr,length(M_n1_list)))
        ng_vec<-cbind(ng_vec,method_vec)
    }
    
    
    
    colnames(BAR_DAT)<-c('Number of detected DE genes','Adjusted P-values threshold','Normalization method')
    
    BAR_DAT<-as.data.frame(BAR_DAT)
    
    BAR_DAT[,1]<-as.numeric(as.character(BAR_DAT[,1]))
    BAR_DAT[,2]<-factor(BAR_DAT[,2],levels=unique(BAR_DAT[,2]))
    BAR_DAT[,3]<-factor(BAR_DAT[,3],levels=unique(BAR_DAT[,3]))
    
    num_gene<-length(M_n1_list[[1]])
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
    cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% methodnames)]
    
    
    N1_BATCH_FIG<-ggplot(data=BAR_DAT, aes(x=BAR_DAT[,2], y=BAR_DAT[,1]/num_gene, fill=BAR_DAT[,3])) +
        geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
        #geom_text(aes(label=round(BAR_DAT[,1]/num_gene,3)), vjust=0, color="black", position = position_dodge(0.9), size=geom_text_size)+
        labs(x = 'Adjusted P values threshold',y='False positive rates',fill='Normalization methods')+ggtitle("") +
        #scale_fill_brewer(palette="Paired")+
        scale_fill_manual(values=cbbPalette[-1] )+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size =textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.position="top",legend.key.size = unit(legend_key_size,"line"))
    return(N1_BATCH_FIG)
}






PCA_FUN<-function(pca_try ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,textsize=4,pointsize=0.5,legendpointsize=1,legend_key_size=1,TITLE=''){

    pr_v<-summary(pca_try )$importance[2,]*100
    percentage <- paste( colnames(pca_try$x), "(", paste( as.character(pr_v), "%", ")", sep="") )
    dat<-pca_try$x[,1:2]
    library(ggplot2)
    BATCH<-ggplot(data=as.data.frame(dat),aes(x=PC1,y=PC2))+geom_point(aes(color=as.factor(LABEL_INDIVIDUAL),shape=as.factor(LABEL_REP)),size=pointsize)+labs(color='Individual',shape='Replicate') + 
        xlab(percentage[1]) + 
        ylab(percentage[2])+
        ggtitle(TITLE)+
        guides(color = guide_legend(override.aes = list(size=legendpointsize)),shape = guide_legend(override.aes = list(size=legendpointsize)))+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size =textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,
              legend.position =c(0.8,0.8),
              #legend.position ='top',
              legend.key.size = unit(legend_key_size,"line"))
    return(BATCH)
}


PCA_FUN2<-function(pca_try ,LABEL_INDIVIDUAL=LABEL_INDIVIDUAL,LABEL_REP=LABEL_REP,TITLE='',legend.x="topleft"){
    
    pr_v<-summary(pca_try )$importance[2,]*100
    percentage <- paste( colnames(pca_try$x), "(", paste( as.character(pr_v), "%", ")", sep="") )
    
test=pca_try$x[,1:2]
indi=unique(LABEL_INDIVIDUAL)
j=c(1,1,1,2,2,2,3,3,3)
repl=c(1,2,3,1,2,3,1,2,3)
plot(test[,1],test[,2],cex=0.1,col="white",xlab='',ylab='')
my.col=c("dark red","red","orangered","dark green","green4","pale green","dark blue","blue","light blue")
for(i in 1:9)
{
    toDo=which((grepl(indi[j[i]],row.names(test))==T)&(LABEL_REP == repl[i]))
    points(test[toDo,1],test[toDo,2],col=my.col[i],pch=22,cex=0.5)
}
legend(x=legend.x,legend=indi,text.col = c("red","green4","blue"),bty="n",cex=0.8)

title(xlab=percentage[1],ylab=percentage[2],cex.main=cex,line=line,cex.lab=cex.lab)
title(TITLE,cex.main=cex,line=line_title,cex.lab=cex.lab)

}
