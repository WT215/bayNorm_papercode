#######smFISH############
#smFISH data was kindly provided by Grun et al.

#1
Notch1_Pou5f1_2i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Notch1A594_Oct4Cy5_Area_2i.csv",header=T)
#2
Pou5f1_Klf4_2i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Oct4_Klf4_cellsize_2i.csv",header=T)
#3
Pou5f1_Klf4_J1S<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Oct4TMR_Klf4Cy5_J1S.csv",header=T)

#4
Pou5f1_Notch1_S<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Oct4TMR_Notch1Cy5_Area_S.csv",header=T)


#5
Sohlh2_Hormad1_2i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Sohlh2A594_Hormad1Cy5_Area_2i.csv",header=T)
#6
Sohlh2_Hormad1_S<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Sohlh2A594_Hormad1Cy5_Area_S.csv",header=T)
#7
Sox2_Pou5f1_2i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Sox2_Oct4_cellsize_2i.csv",header=T)
#8
Sox2_j1s<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Sox2A594_cellarea_j1s.csv",header=T)

#9
Stag3_Gli2_2i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Stag3A594_Gli2Cy5_Area_2i.csv",header=T)
#10
Stag3_Gli2_S<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Stag3A594_Gli2Cy5_Area_S.csv",header=T)

#11
Tpx2_Pcna_J1S<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Tpx2A594_PCNACy5_J1S.csv",header=T)
#12
Tpx2_Pcna_J12i<-read.csv(file="E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/Tpx2A594_PCNACy5_J12i.csv",header=T)

Grun_2i_list<-list(Notch1_Pou5f1_2i=Notch1_Pou5f1_2i,Pou5f1_Klf4_2i=Pou5f1_Klf4_2i,Sohlh2_Hormad1_2i=Sohlh2_Hormad1_2i,Sox2_Pou5f1_2i=Sox2_Pou5f1_2i,Stag3_Gli2_2i=Stag3_Gli2_2i,Tpx2_Pcna_J12i=Tpx2_Pcna_J12i)

# par(mfrow=c(3,4))
# for(i in 1:6){
#   for(j in 1:2){
#     plot(Grun_2i_list[[i]][,3],Grun_2i_list[[i]][,j],xlab='area',ylab='smFISH: gene number',main=paste('2i',colnames(Grun_2i_list[[i]])[j],'cor=',round(cor(Grun_2i_list[[i]][,3],Grun_2i_list[[i]][,j]),4)),pch=16)
#   }
# }

#Grun_2i_namelist<-list(Notch1_Pou5f1_2i=c("Notch1","Pou5f1"),Pou5f1_Klf4_2i=c("Pou5f1","Klf4"),Sohlh2_Hormad1_2i=c("Sohlh2","Hormad1"),Sox2_Pou5f1_2i=c("Sox2","Pou5f1"),Stag3_Gli2_2i=c("Stag3","Gli2"),Tpx2_Pcna_J12i=c("Tpx2","Pcna"))

#Sox2_j1s just one gene:Grun_S_list[4]
Grun_S_list<-list(Pou5f1_Klf4_J1S=Pou5f1_Klf4_J1S,Pou5f1_Notch1_S=Pou5f1_Notch1_S,Sohlh2_Hormad1_S=Sohlh2_Hormad1_S,Sox2_j1s=Sox2_j1s,Stag3_Gli2_S=Stag3_Gli2_S,Tpx2_Pcna_J1S=Tpx2_Pcna_J1S)
names(Grun_S_list)


par(mfrow=c(3,4))
for(i in 1:6){
    if(i!=4){
        
        for(j in 1:2){
            plot(Grun_S_list[[i]][,3],Grun_S_list[[i]][,j],xlab='area',ylab='smFISH: gene number',main=paste('serum',colnames(Grun_S_list[[i]])[j],'cor=',round(cor(Grun_S_list[[i]][,3],Grun_S_list[[i]][,j]),4)),pch=16)
        }
        
    }else if(i==4){
        plot(Grun_S_list[[i]][,2],Grun_S_list[[i]][,1],xlab='area',ylab='smFISH: gene number',main=paste('serum',colnames(Grun_S_list[[i]])[j],'cor=',round(cor(Grun_S_list[[i]][,2],Grun_S_list[[i]][,1]),4)),pch=16)
        
    }
}



interG<-c("Ccna2","Ccnb1","Ccnd1","Gapdh","Pou5f1","Sox2","Pcna",'Klf4',"Sohlh2", "Notch1", "Gli2" ,"Stag3","Tpx2")

length(Grun_2i_list)
length(Grun_S_list)

#####normalized FISH########
Fig3Gene<-c("Sohlh2","Notch1","Gli2","Stag3","Tpx2",'Klf4',"Pcna","Sox2","Pou5f1")
library(foreach)
Grun_2i_norm_list<-list()

smFISH_2i<-foreach(i=1:6,.combine=c)%do%{
    normfa<-Grun_2i_list[[i]][,3]/mean(Grun_2i_list[[i]][,3])
    qq<-apply(Grun_2i_list[[i]],2,function(x){x/normfa})
    Grun_2i_norm_list[[i]]<-qq
    
    cv<-apply(qq,2,function(x){sd(x)/mean(x)})[1:2]
    names(cv)<-colnames(Grun_2i_list[[i]])[c(1,2)]
    return(cv)
}

names(Grun_2i_norm_list)<-names(Grun_2i_list)
summary(Grun_2i_list$Notch1_Pou5f1_2i[,2])



Grun_S_norm_list<-list()
smFISH_S<-foreach(i=1:6,.combine=c)%do%{
    if(i!=4){
        normfa<-Grun_S_list[[i]][,3]/mean(Grun_S_list[[i]][,3])
        qq<-apply(Grun_S_list[[i]],2,function(x){x/normfa})
        Grun_S_norm_list[[i]]<-qq
        cv<-apply(qq,2,function(x){sd(x)/mean(x)})[1:2]
        names(cv)<-colnames(Grun_S_list[[i]])[c(1,2)]
    }else if (i==4){
        normfa<-Grun_S_list[[i]][,2]/mean(Grun_S_list[[i]][,2])
        qq<-apply(Grun_S_list[[i]],2,function(x){x/normfa})
        Grun_S_norm_list[[i]]<-qq
        cv<-apply(qq,2,function(x){sd(x)/mean(x)})[1]
        names(cv)<-colnames(Grun_S_list[[i]])[1]
    }
    return(cv)
}
names(Grun_S_norm_list)<-names(Grun_S_list)



#integrate duplicated genes for 2i
temp2i_mean<-mean(smFISH_2i[which(names(smFISH_2i)=='Pou5f1')])
smFISH_2i_v2<-smFISH_2i[-which(names(smFISH_2i)=='Pou5f1')]
smFISH_2i_v2<-c(smFISH_2i_v2,Pou5f1=temp2i_mean)
smFISH_2i_v2<-smFISH_2i_v2[Fig3Gene]

#integrate duplicated genes for S
tempS_mean<-mean(smFISH_S[which(names(smFISH_S)=='Pou5f1')])
smFISH_S_v2<-smFISH_S[-which(names(smFISH_S)=='Pou5f1')]
smFISH_S_v2<-c(smFISH_S_v2,Pou5f1=tempS_mean)
smFISH_S_v2<-smFISH_S_v2[Fig3Gene]


#only Pou5f1_Klf4_2i used for Pou5f1
LIST_2i<-Grun_2i_norm_list
Fish_2iplot_list<-list(LIST_2i$Sohlh2_Hormad1_2i[,1],LIST_2i$Notch1_Pou5f1_2i[,1],LIST_2i$Stag3_Gli2_2i[,2],LIST_2i$Stag3_Gli2_2i[,1],LIST_2i$Tpx2_Pcna_J12i[,1],LIST_2i$Pou5f1_Klf4_2i[,2],LIST_2i$Tpx2_Pcna_J12i[,2],LIST_2i$Sox2_Pou5f1_2i[,1],LIST_2i$Pou5f1_Klf4_2i[,1])
names(Fish_2iplot_list)<-Fig3Gene

#only Pou5f1_Klf4_J1S used for Pou5f1
LIST_serum<-Grun_S_norm_list
Fish_serumplot_list<-list(LIST_serum$Sohlh2_Hormad1_S[,1],LIST_serum$Pou5f1_Notch1_S[,2],LIST_serum$Stag3_Gli2_S[,2],LIST_serum$Stag3_Gli2_S[,1],LIST_serum$Tpx2_Pcna_J1S[,1],LIST_serum$Pou5f1_Klf4_J1S[,2],LIST_serum$Tpx2_Pcna_J1S[,2],LIST_serum$Sox2_j1s[,1],LIST_serum$Pou5f1_Klf4_J1S[,1])
names(Fish_serumplot_list)<-Fig3Gene

#######end of preparing norm FISH########

#save.image("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/smFISH_norm_load.RData")
load("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/Grun_2014_smFISH/smFISH_norm_load.RData")
