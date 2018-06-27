BULK<-read.table(file="E:/RNAseqProject/Soumillon_2014/GSE53638_D3_Bulk_UMI.dat")
D3_UMI<-read.table("E:/RNAseqProject/Soumillon_2014/GSE53638_D3_UMI.dat")


D3T0<-D3_UMI[,grep(colnames(D3_UMI),pattern='D3T0')]
D3T7<-D3_UMI[,grep(colnames(D3_UMI),pattern='D3T7')]
D3_used<-cbind(D3T0,D3T7)
CONDITION=c(rep(1,dim(D3T0)[2]),rep(2,dim(D3T7)[2]))
names(CONDITION)<-colnames(D3_used)
dim(D3_used)
gr<-c(dim(D3T0)[2],dim(D3T7)[2])
clSum_qun<-quantile(colSums(D3_used),c(0.05,0.95))
D3_used2<-D3_used[,-c(which(colSums(D3_used)<clSum_qun[1] | colSums(D3_used)>clSum_qun[2]))]
D3_used3<-D3_used2[which(rowMeans(D3_used)>0.05),]


#Prepare for the benchmark genes
log2FC<-log2(BULK[rownames(D3_used3),3]+1)-log2(BULK[rownames(D3_used3),1]+1)
names(log2FC)<-rownames(D3_used3)
DE_1000<-names(sort(log2FC,decreasing = T)[1:1000])
DE_TRUE_LABEL<-rep(0,dim(D3_used3)[1])
names(DE_TRUE_LABEL)<-rownames(D3_used3)
DE_TRUE_LABEL[DE_1000]<-1





#subsampling D3_used3 datasets:

#########sample 100 200 400##########
CONDITION_used<-CONDITION[colnames(D3_used3)]
table(CONDITION_used)
D3TO<-CONDITION_used[CONDITION_used==1]
D3T7<-CONDITION_used[CONDITION_used==2]
#####balance group#######
set.seed(12300)
#100
D3T0_100_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),100))
}
D3T7_100_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),100))
}

#200 
D3T0_200_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),200))
}
D3T7_200_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),200))
}

#400 
D3T0_400_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),400))
}
D3T7_400_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),400))
}


save(D3T0_100_gr,D3T7_100_gr,D3T0_200_gr,D3T7_200_gr,D3T0_400_gr,D3T7_400_gr,file="E:/RNAseqProject/Soumillon_2014/D3T_smallsamples.RData")


#########sample 20 50 80##########


CONDITION_used<-CONDITION[colnames(D3_used3)]
table(CONDITION_used)
D3TO<-CONDITION_used[CONDITION_used==1]
D3T7<-CONDITION_used[CONDITION_used==2]
#####balance group#######
set.seed(131256)
#20
D3T0_20_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),20))
}
D3T7_20_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),20))
}

#50
D3T0_50_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),50))
}
D3T7_50_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),50))
}

#80
D3T0_80_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3TO),80))
}
D3T7_80_gr<-foreach(i=1:10)%do%{
    return(sample(names(D3T7),80))
}


save(D3T0_20_gr,D3T7_20_gr,D3T0_50_gr,D3T7_50_gr,D3T0_80_gr,D3T7_80_gr,file="E:/RNAseqProject/Soumillon_2014/D3T_smallsamples_v2.RData")


save.image(file="E:/RNAseqProject/Soumillon_2014/Soumillon_2014.RData")

