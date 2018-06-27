#######DE fun#######
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
                
                
            }
        }
        auc_list[[i]]<-auc_vec
        
        
        
    }
    names(auc_list)<-names(LIST) 
    return(auc_list)
    
}


#####load raw data#########
load("E:/RNAseqProject/Soumillon_2014/D3T_smallsamples_v2.RData")
load("Soumillon_norms.RData")


log2FC<-log2(BULK[rownames(D3_used3),3]+1)-log2(BULK[rownames(D3_used3),1]+1)
names(log2FC)<-rownames(D3_used3)
DE_1000<-names(sort(log2FC,decreasing = T)[1:1000])

DE_TRUE_LABEL<-rep(0,dim(D3_used3)[1])
names(DE_TRUE_LABEL)<-rownames(D3_used3)
DE_TRUE_LABEL[DE_1000]<-1


TRUE_LABEL_input<-DE_TRUE_LABEL

TRUE_LABEL_input[TRUE_LABEL_input==0]=3
TRUE_LABEL_input[TRUE_LABEL_input==1]=0
TRUE_LABEL_input[TRUE_LABEL_input==3]=1

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

###bayNorm###########
##20 vs 20
library(bayNorm)
library(abind)
library(foreach)
temp_bay<-abind(aBay_out$Bay_array_list$`Group 1`,aBay_out$Bay_array_list$`Group 2`,along=2)
dim(temp_bay)
bayNorm_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(20,20))
    return(qq)
}



##50 vs 50
bayNorm_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(50,50))
    return(qq)
}



##80 vs 80
bayNorm_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(80,80))
    return(qq)
}


##20 vs 50
bayNorm_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(20,50))
    return(qq)
}



##50 vs 20
bayNorm_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(50,20))
    return(qq)
}



## 20 vs 80
bayNorm_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(20,80))
    return(qq)
}


##80 vs 20
bayNorm_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(80,20))
    return(qq)
}


##50 vs 80
bayNorm_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(50,80))
    return(qq)
}


##80 vs 50
bayNorm_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(80,50))
    return(qq)
}




LIST_bayNorm<-list(bayNorm_20_20=bayNorm_20_20,bayNorm_50_50=bayNorm_50_50,bayNorm_80_80=bayNorm_80_80,bayNorm_20_50=bayNorm_20_50,bayNorm_20_80=bayNorm_20_80,bayNorm_50_20=bayNorm_50_20,bayNorm_80_20=bayNorm_80_20,bayNorm_50_80=bayNorm_50_80,bayNorm_80_50=bayNorm_80_50)

save(LIST_bayNorm,file="E:/RNAseqProject/Soumillon_2014/smallgroup_bayNorm_v2.RData")
load("E:/RNAseqProject/Soumillon_2014/smallgroup_bayNorm_v2.RData")

DE_bay<-DE_smallfun(LIST=LIST_bayNorm,TRUELABEL=TRUE_LABEL_input)
DE_bay
unlist(lapply(DE_bay,mean))
#########saver############
#load("E:/RNAseqProject/Soumillon_2014/Soumillon_SAVER.RData")
#load("E:/RNAseqProject/Soumillon_2014/saver_beta/Soumillon_saver_BETA.RData")
load("E:/RNAseqProject/Soumillon_2014/saver_default/Soumillon_saver_default.RData")
##100 vs 100
library(bayNorm)
library(abind)
temp_SAVER<-abind(saver_gr1_a,saver_gr2_a,along=2)
##20 vs 20
library(foreach)
SAVER_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(20,20))
    cat("\014")
    return(qq)
}

cat("\014")

##50 vs 50
SAVER_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(50,50))
    cat("\014")
    return(qq)
}

cat("\014")

##80 vs 80
SAVER_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(80,80))
    cat("\014")
    return(qq)
}
cat("\014")

##20 vs 50
SAVER_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(20,50))
    cat("\014")
    return(qq)
}
cat("\014")


##50 vs 20
SAVER_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(50,20))
    cat("\014")
    return(qq)
}
cat("\014")


## 20 vs 80
SAVER_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(20,80))
    cat("\014")
    return(qq)
}

cat("\014")
##80 vs 20
SAVER_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(80,20))
    cat("\014")
    return(qq)
}
cat("\014")

##50 vs 80
SAVER_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(50,80))
    cat("\014")
    return(qq)
}

cat("\014")
##80 vs 50
SAVER_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(80,50))
    cat("\014")
    return(qq)
}

cat("\014")


LIST_SAVER<-list(SAVER_20_20=SAVER_20_20,SAVER_50_50=SAVER_50_50,SAVER_80_80=SAVER_80_80,SAVER_20_50=SAVER_20_50,SAVER_20_80=SAVER_20_80,SAVER_50_20=SAVER_50_20,SAVER_80_20=SAVER_80_20,SAVER_50_80=SAVER_50_80,SAVER_80_50=SAVER_80_50)

DE_SAVER<-DE_smallfun(LIST=LIST_SAVER,TRUELABEL=TRUE_LABEL_input)
DE_SAVER
unlist(lapply(DE_SAVER,mean))

save(LIST_SAVER,file="E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER_default_v2.RData")

load("E:/RNAseqProject/Soumillon_2014/smallgroup_SAVER_v2.RData")

#####SCnorm##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##20 vs 20
library(bayNorm)
library(abind)


SCnorm_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(20,20))
    return(qq)
}
cat("\014")  


##50 vs 50


SCnorm_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(50,50))
    return(qq)
}
cat("\014")  


##80 vs 80

SCnorm_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(80,80))
    return(qq)
}
cat("\014")  

##20 vs 50

SCnorm_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(20,50))
    return(qq)
}

cat("\014")  

##50 vs 20


SCnorm_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(50,20))
    return(qq)
}
cat("\014")  


## 20vs 80

SCnorm_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(20,80))
    return(qq)
}
cat("\014")  

##80 vs 20


SCnorm_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(80,20))
    return(qq)
}
cat("\014")  

##50 vs 80

SCnorm_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(50,80))
    return(qq)
}
cat("\014")  
##80 vs 50


SCnorm_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(80,50))
    return(qq)
}



LIST_SCnorm<-list(SCnorm_20_20=SCnorm_20_20,SCnorm_50_50=SCnorm_50_50,SCnorm_80_80=SCnorm_80_80,SCnorm_20_50=SCnorm_20_50,SCnorm_20_80=SCnorm_20_80,SCnorm_50_20=SCnorm_50_20,SCnorm_80_20=SCnorm_80_20,SCnorm_50_80=SCnorm_50_80,SCnorm_80_50=SCnorm_80_50)


save(LIST_SCnorm,file="E:/RNAseqProject/Soumillon_2014/smallgroup_SCnorm_v2.RData")

#####Scaling##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")
log2FC<-log2(BULK[rownames(D3_used3),3]+1)-log2(BULK[rownames(D3_used3),1]+1)
names(log2FC)<-rownames(D3_used3)
DE_1000<-names(sort(log2FC,decreasing = T)[1:1000])
source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##20 vs 20
library(bayNorm)
library(abind)


RB_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(20,20))
    return(qq)
}

cat("\014")  

##50 vs 50
RB_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(50,50))
    return(qq)
}


cat("\014")  
##80 vs 80
RB_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(80,80))
    return(qq)
}

cat("\014")  
##20 vs 50
RB_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(20,50))
    return(qq)
}
cat("\014")  


##50 vs 20
RB_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(50,20))
    return(qq)
}

cat("\014")  

## 20vs 80
RB_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(20,80))
    return(qq)
}

cat("\014")  
##80 vs 20
RB_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(80,20))
    return(qq)
}
cat("\014")  

##50 vs 80
RB_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(50,80))
    return(qq)
}
cat("\014")  

##80 vs 50
RB_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(80,50))
    return(qq)
}






LIST_Scaling<-list(RB_20_20=RB_20_20,RB_50_50=RB_50_50,RB_80_80=RB_80_80,RB_20_50=RB_20_50,RB_20_80=RB_20_80,RB_50_20=RB_50_20,RB_80_20=RB_80_20,RB_50_80=RB_50_80,RB_80_50=RB_80_50)
save(LIST_Scaling,file="E:/RNAseqProject/Soumillon_2014/smallgroup_Scaling_v2.RData")




#####scImpute##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")


source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##20 vs 20
library(bayNorm)
library(abind)


scImpute_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(20,20))
    return(qq)
}

cat("\014")  

##50 vs 50
scImpute_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(50,50))
    return(qq)
}
cat("\014")  


##80 vs 80
scImpute_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(80,80))
    return(qq)
}
cat("\014")  

##20 vs 50
scImpute_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(20,50))
    return(qq)
}
cat("\014")  


##50 vs 20
scImpute_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(50,20))
    return(qq)
}
cat("\014")  


## 20vs 80
scImpute_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(20,80))
    return(qq)
}
cat("\014")  

##80 vs 20
scImpute_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(80,20))
    return(qq)
}
cat("\014")  

##50 vs 80
scImpute_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(50,80))
    return(qq)
}
cat("\014")  

##80 vs 50
scImpute_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(80,50))
    return(qq)
}


LIST_scImpute<-list(scImpute_20_20=scImpute_20_20,scImpute_50_50=scImpute_50_50,scImpute_80_80=scImpute_80_80,scImpute_20_50=scImpute_20_50,scImpute_20_80=scImpute_20_80,scImpute_50_20=scImpute_50_20,scImpute_80_20=scImpute_80_20,scImpute_50_80=scImpute_50_80,scImpute_80_50=scImpute_80_50)
save(LIST_scImpute,file="E:/RNAseqProject/Soumillon_2014/smallgroup_scImpute_v2.RData")






#####MAGIC##############
#load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")
##20 vs 20
library(bayNorm)
library(abind)


MAGIC_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(20,20))
    return(qq)
}



##50 vs 50
MAGIC_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(50,50))
    return(qq)
}



##80 vs 80
MAGIC_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(80,80))
    return(qq)
}


##20 vs 50
MAGIC_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(20,50))
    return(qq)
}



##50 vs 20
MAGIC_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(50,20))
    return(qq)
}



## 20vs 80
MAGIC_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(20,80))
    return(qq)
}


##80 vs 20
MAGIC_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(80,20))
    return(qq)
}


##50 vs 80
MAGIC_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(50,80))
    return(qq)
}


##80 vs 50
MAGIC_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(80,50))
    return(qq)
}


LIST_MAGIC<-list(MAGIC_20_20=MAGIC_20_20,MAGIC_50_50=MAGIC_50_50,MAGIC_80_80=MAGIC_80_80,MAGIC_20_50=MAGIC_20_50,MAGIC_20_80=MAGIC_20_80,MAGIC_50_20=MAGIC_50_20,MAGIC_80_20=MAGIC_80_20,MAGIC_50_80=MAGIC_50_80,MAGIC_80_50=MAGIC_80_50)
save(LIST_MAGIC,file="E:/RNAseqProject/Soumillon_2014/smallgroup_MAGIC_v2.RData")



#####DCA##############
#load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")
##20 vs 20
library(bayNorm)
library(abind)


DCA_20_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(20,20))
    return(qq)
}



##50 vs 50
DCA_50_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(50,50))
    return(qq)
}



##80 vs 80
DCA_80_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(80,80))
    return(qq)
}


##20 vs 50
DCA_20_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(20,50))
    return(qq)
}



##50 vs 20
DCA_50_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(50,20))
    return(qq)
}



## 20vs 80
DCA_20_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_20_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(20,80))
    return(qq)
}


##80 vs 20
DCA_80_20<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_20_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(80,20))
    return(qq)
}


##50 vs 80
DCA_50_80<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_50_gr[[i]],D3T7_80_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(50,80))
    return(qq)
}


##80 vs 50
DCA_80_50<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_80_gr[[i]],D3T7_50_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(80,50))
    return(qq)
}


LIST_DCA<-list(DCA_20_20=DCA_20_20,DCA_50_50=DCA_50_50,DCA_80_80=DCA_80_80,DCA_20_50=DCA_20_50,DCA_20_80=DCA_20_80,DCA_50_20=DCA_50_20,DCA_80_20=DCA_80_20,DCA_50_80=DCA_50_80,DCA_80_50=DCA_80_50)
save(LIST_DCA,file="E:/RNAseqProject/Soumillon_2014/smallgroup_DCA_v2.RData")
