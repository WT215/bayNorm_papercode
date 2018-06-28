CV_fun<-function(data){
    
    if(length(dim(data))!=3){
        data_sd<-apply(data,1,sd)
        cv<-data_sd  /rowMeans(data)
    } else if(length(dim(data))==3){
        MuTest<-apply(data,1,mean)
        SdTest<-apply(apply(data, c(1,3), sd), 1, mean)
        cv<-SdTest/MuTest
    }
    
    return(cv)
}

SIM_ana_fun<-function(DATA_LIST){
    library(foreach)
    library(doSNOW)
    library(parallel)
    library(ineq)
    cluster = makeCluster(5, type = "SOCK")
    registerDoSNOW(cluster)
    getDoParWorkers()
    
    iterations <-length(DATA_LIST)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    
    qq<-foreach(i=1:length(DATA_LIST),.packages = 'ineq',.export = 'CV_fun',.options.snow = opts)%dopar%{
        cv_out<-CV_fun(DATA_LIST[[i]])
        mean_out<-rowMeans(DATA_LIST[[i]])
        
        gini_out<- apply(DATA_LIST[[i]],1,Gini)
        
        re_li<-list(mean=mean_out,cv=cv_out,gini=gini_out)
        return(re_li)
    }
    close(pb)
    stopCluster(cluster)
    names(qq)<-names(DATA_LIST)
    
    return(qq)
    
}

mean_data_fun<-function(temp_list,TRUE_MU,textsize=10,lwd=1,outlier.size=0.005){
    

mean_dat<-foreach(i=1:length(temp_list),.combine=rbind)%:%
    foreach(j=1:length(temp_list[[i]]),.combine=rbind)%do%{
        ll<-length(temp_list[[i]][[j]]$mean)
        if(i==1){
            qtemp<-log2(temp_list[[i]][[j]]$mean/TRUE_MU)
        }else{
            qtemp<-log2(temp_list[[i]][[j]]$mean/TRUE_MU)
        }
        
        
        qq<-cbind(qtemp,rep(names(temp_list[[i]])[j],ll),rep(paste('Group',i),ll))  
    }
colnames(mean_dat)<-c("log2 ratio of mean: scRNA-seq / True count",'Method','Group')

mean_dat<-as.data.frame(mean_dat)
mean_dat[,1]<-as.numeric(as.character(mean_dat[,1]))
mean_dat[,2]<-factor(mean_dat[,2],levels=unique(mean_dat[,2]))
mean_dat$Group<-factor(mean_dat$Group,levels=unique(mean_dat$Group))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


MEAN_PLOT=ggplot(mean_dat,aes(x=Method,y=mean_dat[,1],fill=Group))+
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9,binwidth=0.15)+
    scale_fill_manual(values=cbbPalette[-1])+
    labs(x='',y="log2 ratio of mean: scRNA-seq / True count")+
    geom_boxplot(outlier.size = outlier.size,outlier.shape = 46,lwd=lwd,outlier.colour =NA)+
    geom_hline(yintercept=0,alpha=0.25)+
    geom_hline(yintercept=c(-1,1),lty=2,alpha=0.5)+
    
    ylim(-2.5,2.5)+
    
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
return(list(MEAN_PLOT=MEAN_PLOT,mean_dat=mean_dat))
}




cv_data_fun<-function(temp_list,TRUE_CV_g1,TRUE_CV_g2,textsize=10,lwd=1,outlier.size=0.005){
    
    
    cv_dat<-foreach(i=1:length(temp_list),.combine=rbind)%:%
        foreach(j=1:length(temp_list[[i]]),.combine=rbind)%do%{
            ll<-length(temp_list[[i]][[j]]$mean)
            if(i==1){
                qtemp<-log2(temp_list[[i]][[j]]$cv/TRUE_CV_g1)
            }else{
                qtemp<-log2(temp_list[[i]][[j]]$cv/TRUE_CV_g2)
            }
            
            
            qq<-cbind(qtemp,rep(names(temp_list[[i]])[j],ll),rep(paste('Group',i),ll))  
        }
    colnames(cv_dat)<-c("log2 ratio of cv: scRNA-seq / True count",'Method','Group')
    
    cv_dat<-as.data.frame(cv_dat)
    cv_dat[,1]<-as.numeric(as.character(cv_dat[,1]))
    cv_dat[,2]<-factor(cv_dat[,2],levels=unique(cv_dat[,2]))
    cv_dat$Group<-factor(cv_dat$Group,levels=unique(cv_dat$Group))
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    
    CV_PLOT=ggplot(cv_dat,aes(x=Method,y=cv_dat[,1],fill=Group))+
        scale_fill_manual(values=cbbPalette[-1])+
        labs(x='',y="log2 ratio of cv: scRNA-seq / True count")+
        geom_boxplot(outlier.size =outlier.size,outlier.shape = 46,lwd=lwd,outlier.colour =NA)+
        geom_hline(yintercept=0,alpha=0.25)+
        geom_hline(yintercept=c(-1,1),lty=2,alpha=0.5)+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
    return(list(CV_PLOT=CV_PLOT,cv_dat=cv_dat))
}



gini_data_fun<-function(temp_list,TRUE_Gini_g1,TRUE_Gini_g2,textsize=10,lwd=1,outlier.size=0.005){
    
    
    gini_dat<-foreach(i=1:length(temp_list),.combine=rbind)%:%
        foreach(j=1:length(temp_list[[i]]),.combine=rbind)%do%{
            ll<-length(temp_list[[i]][[j]]$mean)
            if(i==1){
                qtemp<-log2(temp_list[[i]][[j]]$gini/TRUE_Gini_g1)
            }else{
                qtemp<-log2(temp_list[[i]][[j]]$gini/TRUE_Gini_g2)
            }
            
            
            qq<-cbind(qtemp,rep(names(temp_list[[i]])[j],ll),rep(paste('Group',i),ll))  
        }
    colnames(gini_dat)<-c("log2 ratio of Gini: scRNA-seq / True count",'Method','Group')
    
    gini_dat<-as.data.frame(gini_dat)
    gini_dat[,1]<-as.numeric(as.character(gini_dat[,1]))
    gini_dat[,2]<-factor(gini_dat[,2],levels=unique(gini_dat[,2]))
    gini_dat$Group<-factor(gini_dat$Group,levels=unique(gini_dat$Group))
    
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    
    Gini_PLOT=ggplot(gini_dat,aes(x=Method,y=gini_dat[,1],fill=Group))+
        scale_fill_manual(values=cbbPalette[-1])+
        labs(x='',y="log2 ratio of Gini: scRNA-seq / True count")+
        geom_boxplot(outlier.size = outlier.size,outlier.shape = 46,lwd=lwd,outlier.colour =NA)+
        geom_hline(yintercept=0,alpha=0.25)+
        geom_hline(yintercept=c(-1,1),lty=2,alpha=0.5)+
        theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(0.8,"line"),legend.position ='top')
    return(list(Gini_PLOT=Gini_PLOT,gini_dat=gini_dat))
}