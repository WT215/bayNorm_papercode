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

load("E:/RNAseqProject/Soumillon_2014/D3T_smallsamples.RData")
load("Soumillon_2014.RData")
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


source("E:/RNAseqProject/MANY_DE_FUN.R")


###bayNorm###########
##100 vs 100
library(bayNorm)
library(abind)
library(foreach)
temp_bay<-abind(aBay_out$Bay_array_list$`Group 1`,aBay_out$Bay_array_list$`Group 2`,along=2)
dim(temp_bay)
bayNorm_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(100,100))
    return(qq)
}



##200 vs 200
bayNorm_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(200,200))
    return(qq)
}



##400 vs 400
bayNorm_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(400,400))
    return(qq)
}


##100 vs 200
bayNorm_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(100,200))
    return(qq)
}



##200 vs 100
bayNorm_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(200,100))
    return(qq)
}



## 100 vs 400
bayNorm_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(100,400))
    return(qq)
}


##400 vs 100
bayNorm_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
bayNorm_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(200,400))
    return(qq)
}


##400 vs 200
bayNorm_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_bay[,cellused,],NumCells=c(400,200))
    return(qq)
}




LIST_bayNorm<-list(bayNorm_100_100=bayNorm_100_100,bayNorm_200_200=bayNorm_200_200,bayNorm_400_400=bayNorm_400_400,bayNorm_100_200=bayNorm_100_200,bayNorm_100_400=bayNorm_100_400,bayNorm_200_100=bayNorm_200_100,bayNorm_400_100=bayNorm_400_100,bayNorm_200_400=bayNorm_200_400,bayNorm_400_200=bayNorm_400_200)

save(LIST_bayNorm,file="E:/RNAseqProject/Soumillon_2014/smallgroup_bayNorm.RData")

DE_bay<-DE_smallfun(LIST=LIST_bayNorm,TRUELABEL=TRUE_LABEL_input)
DE_bay
unlist(lapply(DE_bay,mean))
#########saver############
#load("E:/RNAseqProject/Soumillon_2014/saver_beta/Soumillon_saver_BETA.RData")
load("E:/RNAseqProject/Soumillon_2014/saver_default/Soumillon_saver_default.RData")


length(which(rowSums(saver_gr2$estimate)==0))
length(which(cbind(saver_gr1$estimate,saver_gr2$estimate)==0))
length(which(D3_used3==0))
length(which(rowSums(saver_gr1_a)==0))
dim(saver_gr1_a)
##100 vs 100
library(bayNorm)
library(abind)
temp_SAVER<-abind(saver_gr1_a,saver_gr2_a,along=2)
dim(temp_bay)
SAVER_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(100,100))
    return(qq)
}

cat("\014")  

##200 vs 200
SAVER_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(200,200))
    return(qq)
}

cat("\014")  

##400 vs 400
SAVER_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(400,400))
    return(qq)
}

cat("\014")  
##100 vs 200
SAVER_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(100,200))
    return(qq)
}

cat("\014")  

##200 vs 100
SAVER_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(200,100))
    return(qq)
}

cat("\014")  

## 100 vs 400
SAVER_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(100,400))
    return(qq)
}
cat("\014")  

##400 vs 100
SAVER_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(400,100))
    return(qq)
}
cat("\014")  

##200 vs 400
SAVER_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(200,400))
    return(qq)
}
cat("\014")  

##400 vs 200
SAVER_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=temp_SAVER[,cellused,],NumCells=c(400,200))
    return(qq)
}
cat("\014")  

LIST_SAVER<-list(SAVER_100_100=SAVER_100_100,SAVER_200_200=SAVER_200_200,SAVER_400_400=SAVER_400_400,SAVER_100_200=SAVER_100_200,SAVER_100_400=SAVER_100_400,SAVER_200_100=SAVER_200_100,SAVER_400_100=SAVER_400_100,SAVER_200_400=SAVER_200_400,SAVER_400_200=SAVER_400_200)

save(LIST_SAVER,file="E:/RNAseqProject/Soumillon_2014/saver_default/smallgroup_SAVER_default.RData")



#####SCnorm##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##100 vs 100
library(bayNorm)
library(abind)

SCnorm_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(100,100))
    return(qq)
}



##200 vs 200
SCnorm_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(200,200))
    return(qq)
}



##400 vs 400
SCnorm_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(400,400))
    return(qq)
}


##100 vs 200
SCnorm_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(100,200))
    return(qq)
}



##200 vs 100
SCnorm_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(200,100))
    return(qq)
}



## 100 vs 400
SCnorm_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(100,400))
    return(qq)
}


##400 vs 100
SCnorm_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
SCnorm_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(200,400))
    return(qq)
}


##400 vs 200
SCnorm_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=SCnorm_out@metadata$NormalizedData[,cellused],NumCells=c(400,200))
    return(qq)
}
LIST_SCnorm<-list(SCnorm_100_100=SCnorm_100_100,SCnorm_200_200=SCnorm_200_200,SCnorm_400_400=SCnorm_400_400,SCnorm_100_200=SCnorm_100_200,SCnorm_100_400=SCnorm_100_400,SCnorm_200_100=SCnorm_200_100,SCnorm_400_100=SCnorm_400_100,SCnorm_200_400=SCnorm_200_400,SCnorm_400_200=SCnorm_400_200)
save(LIST_SCnorm,file="E:/RNAseqProject/Soumillon_2014/smallgroup_SCnorm.RData")

#####Scaling##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")
log2FC<-log2(BULK[rownames(D3_used3),3]+1)-log2(BULK[rownames(D3_used3),1]+1)
names(log2FC)<-rownames(D3_used3)
DE_1000<-names(sort(log2FC,decreasing = T)[1:1000])
source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##100 vs 100
library(bayNorm)
library(abind)

RB_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(100,100))
    return(qq)
}



##200 vs 200
RB_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(200,200))
    return(qq)
}



##400 vs 400
RB_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(400,400))
    return(qq)
}


##100 vs 200
RB_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(100,200))
    return(qq)
}



##200 vs 100
RB_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(200,100))
    return(qq)
}



## 100 vs 400
RB_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(100,400))
    return(qq)
}


##400 vs 100
RB_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
RB_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(200,400))
    return(qq)
}


##400 vs 200
RB_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=RB[,cellused],NumCells=c(400,200))
    return(qq)
}
LIST_Scaling<-list(RB_100_100=RB_100_100,RB_200_200=RB_200_200,RB_400_400=RB_400_400,RB_100_200=RB_100_200,RB_100_400=RB_100_400,RB_200_100=RB_200_100,RB_400_100=RB_400_100,RB_200_400=RB_200_400,RB_400_200=RB_400_200)
save(LIST_Scaling,file="E:/RNAseqProject/Soumillon_2014/smallgroup_Scaling.RData")




#####scImpute##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##100 vs 100
library(bayNorm)
library(abind)

scImpute_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(100,100))
    return(qq)
}



##200 vs 200
scImpute_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(200,200))
    return(qq)
}



##400 vs 400
scImpute_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(400,400))
    return(qq)
}


##100 vs 200
scImpute_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(100,200))
    return(qq)
}



##200 vs 100
scImpute_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(200,100))
    return(qq)
}



## 100 vs 400
scImpute_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(100,400))
    return(qq)
}


##400 vs 100
scImpute_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
scImpute_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(200,400))
    return(qq)
}


##400 vs 200
scImpute_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=scImpute_out[,cellused],NumCells=c(400,200))
    return(qq)
}
LIST_scImpute<-list(scImpute_100_100=scImpute_100_100,scImpute_200_200=scImpute_200_200,scImpute_400_400=scImpute_400_400,scImpute_100_200=scImpute_100_200,scImpute_100_400=scImpute_100_400,scImpute_200_100=scImpute_200_100,scImpute_400_100=scImpute_400_100,scImpute_200_400=scImpute_200_400,scImpute_400_200=scImpute_400_200)
save(LIST_scImpute,file="E:/RNAseqProject/Soumillon_2014/smallgroup_scImpute.RData")





#####MAGIC##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")

##100 vs 100
library(bayNorm)
library(abind)

MAGIC_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(100,100))
    return(qq)
}



##200 vs 200
MAGIC_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(200,200))
    return(qq)
}



##400 vs 400
MAGIC_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(400,400))
    return(qq)
}


##100 vs 200
MAGIC_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(100,200))
    return(qq)
}



##200 vs 100
MAGIC_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(200,100))
    return(qq)
}



## 100 vs 400
MAGIC_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(100,400))
    return(qq)
}


##400 vs 100
MAGIC_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
MAGIC_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(200,400))
    return(qq)
}


##400 vs 200
MAGIC_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=MAGIC_out[,cellused],NumCells=c(400,200))
    return(qq)
}
LIST_MAGIC<-list(MAGIC_100_100=MAGIC_100_100,MAGIC_200_200=MAGIC_200_200,MAGIC_400_400=MAGIC_400_400,MAGIC_100_200=MAGIC_100_200,MAGIC_100_400=MAGIC_100_400,MAGIC_200_100=MAGIC_200_100,MAGIC_400_100=MAGIC_400_100,MAGIC_200_400=MAGIC_200_400,MAGIC_400_200=MAGIC_400_200)
save(LIST_MAGIC,file="E:/RNAseqProject/Soumillon_2014/smallgroup_MAGIC.RData")







#####DCA##############
load("E:/RNAseqProject/Soumillon_2014/Soumillon_analysis.RData")
load("E:/RNAseqProject/Soumillon_2014/D3T_smallsamples.RData")
load("E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_Soumillon/DCA_Soumillon.RData")

source("E:/RNAseqProject/MANY_SAVE_PATH.r")
source("E:/RNAseqProject/TSTAT_140817.r")
source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/MANY_NORM_FUN.R")


DCA_out<-DCA_Soumillon
##100 vs 100
library(bayNorm)
library(abind)

DCA_100_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(100,100))
    return(qq)
}

cat("\014")

##200 vs 200
DCA_200_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(200,200))
    return(qq)
}
cat("\014")


##400 vs 400
DCA_400_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(400,400))
    return(qq)
}
cat("\014")

##100 vs 200
DCA_100_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(100,200))
    return(qq)
}


cat("\014")
##200 vs 100
DCA_200_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(200,100))
    return(qq)
}

cat("\014")

## 100 vs 400
DCA_100_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_100_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(100,400))
    return(qq)
}
cat("\014")

##400 vs 100
DCA_400_100<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_100_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(400,100))
    return(qq)
}


##200 vs 400
DCA_200_400<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_200_gr[[i]],D3T7_400_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(200,400))
    return(qq)
}

cat("\014")
##400 vs 200
DCA_400_200<-foreach(i=1:10)%do%{
    print(i)
    cellused<-c(D3T0_400_gr[[i]],D3T7_200_gr[[i]])
    qq<-SCnorm_runMAST3(Data=DCA_out[,cellused],NumCells=c(400,200))
    return(qq)
}
cat("\014")
LIST_DCA<-list(DCA_100_100=DCA_100_100,DCA_200_200=DCA_200_200,DCA_400_400=DCA_400_400,DCA_100_200=DCA_100_200,DCA_100_400=DCA_100_400,DCA_200_100=DCA_200_100,DCA_400_100=DCA_400_100,DCA_200_400=DCA_200_400,DCA_400_200=DCA_400_200)
save(LIST_DCA,file="E:/RNAseqProject/Soumillon_2014/smallgroup_DCA.RData")
