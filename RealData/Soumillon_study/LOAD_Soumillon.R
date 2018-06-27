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

save.image(file="E:/RNAseqProject/Soumillon_2014/Soumillon_2014.RData")

