library(readr)
H1_p24<- read_csv("E:/RNAseqProject/NEWPROJECT_PAPERS/RNAseq_Rpackage_relatedPaper/scRNA_package/SCnorm/SCnorm_dat/H1b5s_EC_p24_93cells.csv")

H1_p24<-as.data.frame(H1_p24)
rownames(H1_p24)<-H1_p24[,1]
H1_p24<-H1_p24[,-1]
dim(H1_p24)


H1_p96<- read_csv("E:/RNAseqProject/NEWPROJECT_PAPERS/RNAseq_Rpackage_relatedPaper/scRNA_package/SCnorm/SCnorm_dat/H1b5s_EC_p96_93cells.csv")
H1_p96<-as.data.frame(H1_p96)
rownames(H1_p96)<-H1_p96[,1]
H1_p96<-H1_p96[,-1]
dim(H1_p96)




########H9######
library(readr)
H9_p24<- read_csv("E:/RNAseqProject/NEWPROJECT_PAPERS/RNAseq_Rpackage_relatedPaper/scRNA_package/SCnorm/SCnorm_dat/H9b3s_EC_p24_91cells.csv")
H9_p24<-as.data.frame(H9_p24)
rownames(H9_p24)<-H9_p24[,1]
H9_p24<-H9_p24[,-1]
dim(H9_p24)


library(readr)
H9_p96<- read_csv("E:/RNAseqProject/NEWPROJECT_PAPERS/RNAseq_Rpackage_relatedPaper/scRNA_package/SCnorm/SCnorm_dat/H9b3s_EC_p96_91cells.csv")
H9_p96<-as.data.frame(H9_p96)
rownames(H9_p96)<-H9_p96[,1]
H9_p96<-H9_p96[,-1]
dim(H9_p24)


#####subset of genes for DE detection########


#Genenames<-rownames(H9_p24)
#Genenames[5292:5383]
ERCC_ind<-seq(5292,5383)

ERCC_H1_p24<-H1_p24[ERCC_ind,]
ERCC_H1_p96<-H1_p96[ERCC_ind,]
ERCC_H9_p24<-H9_p24[ERCC_ind,]
ERCC_H9_p96<-H9_p96[ERCC_ind,]

# ERCC_beta_H1_p24<-colSums(ERCC_H1_p24)/total_ercc_molecules
# ERCC_beta_H1_p96<-colSums(ERCC_H1_p96)/total_ercc_molecules
# ERCC_beta_H9_p24<-colSums(ERCC_H9_p24)/total_ercc_molecules
# ERCC_beta_H9_p96<-colSums(ERCC_H9_p96)/total_ercc_molecules


H1_p24<-H1_p24[-ERCC_ind,]
H1_p96<-H1_p96[-ERCC_ind,]
H9_p24<-H9_p24[-ERCC_ind,]
H9_p96<-H9_p96[-ERCC_ind,]


drop1<-which(rowSums(H1_p24)==0)
drop2<-which(rowSums(H1_p96)==0)
drop3<-which(rowSums(H9_p24)==0)
drop4<-which(rowSums(H9_p96)==0)
drop<-Reduce(intersect, list(drop1,drop2,drop3,drop4))

H1_p24<-H1_p24[-drop,]
H1_p96<-H1_p96[-drop,]
H9_p24<-H9_p24[-drop,]
H9_p96<-H9_p96[-drop,]

# H1_p24<-round(H1_p24)
# H1_p96<-round(H1_p96)
# H9_p24<-round(H9_p24)
# H9_p96<-round(H9_p96)

NZeros_H1_p24 <- apply(H1_p24, 1, function(x) sum(x!=0))
NZeros_H1_p96 <- apply(H1_p96, 1, function(x) sum(x!=0))
# Which genes have at least 10 non-zeros.
NZeros_H1_p24 <- NZeros_H1_p24[which(NZeros_H1_p24 >= 10)]
NZeros_H1_p96 <- NZeros_H1_p96[which(NZeros_H1_p96 >= 10)]
# Consider genes with at least 10 non-zeros.
whichg_H1 <- intersect(names(NZeros_H1_p24), names(NZeros_H1_p96));length(whichg_H1)


NZeros_H9_p24 <- apply(H9_p24, 1, function(x) sum(x!=0))
NZeros_H9_p96 <- apply(H9_p96, 1, function(x) sum(x!=0))
# Which genes have at least 10 non-zeros.
NZeros_H9_p24 <- NZeros_H9_p24[which(NZeros_H9_p24 >= 10)]
NZeros_H9_p96 <- NZeros_H9_p96[which(NZeros_H9_p96 >= 10)]
# Consider genes with at least 10 non-zeros.
whichg_H9 <- intersect(names(NZeros_H9_p24), names(NZeros_H9_p96))
length(whichg_H9)
length(intersect(whichg_H9,whichg_H1))


save.image("RAW_INITIATE.RData")
