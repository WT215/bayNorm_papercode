mydata = read.table("E:/RNAseqProject/tung2017batch/molecules-filter.txt")


DECENT_UP = read.table("E:/RNAseqProject/tung2017batch/DESCENT_DE/tung_bulk_191vs192_up_voom.txt")
DECENT_UP $V1<-as.character(DECENT_UP $V1)
DECENT_DN = read.table("E:/RNAseqProject/tung2017batch/DESCENT_DE/tung_bulk_191vs192_dn_voom.txt")
DECENT_DN $V1<-as.character(DECENT_DN$V1)

length(intersect(DECENT_UP$V1,rownames(mydata)))


dim(N1_1_DAT)

DE_TRUE_LABEL<-rep(0,13058)
names(DE_TRUE_LABEL)<-rownames(mydata)[!ERCC_ind]
dim(mydata)
length(intersect(DECENT_UP$V1,rownames(mydata)))
length(DECENT_UP$V1)
length(intersect(DECENT_DN$V1,rownames(mydata)))
length(DECENT_DN$V1)


DE_TRUE_LABEL[intersect(DECENT_UP$V1,rownames(N1_DAT))]=1
DE_TRUE_LABEL[intersect(DECENT_DN$V1,rownames(N1_DAT))]=1

sum(DE_TRUE_LABEL)


N1_1<-grepl("NA19098.r1", colnames(mydata))
N1_2<-grepl("NA19098.r2", colnames(mydata))
N1_3<-grepl("NA19098.r3", colnames(mydata))

N2_1<-grepl("NA19101.r1", colnames(mydata))
N2_2<-grepl("NA19101.r2", colnames(mydata))
N2_3<-grepl("NA19101.r3", colnames(mydata))

N3_1<-grepl("NA19239.r1", colnames(mydata))
N3_2<-grepl("NA19239.r2", colnames(mydata))
N3_3<-grepl("NA19239.r3", colnames(mydata))

#ERCC
ERCC_ind<-grepl("ERCC", rownames(mydata))
sum(ERCC_ind)
rownames(mydata)[ERCC_ind]





N1_1_DAT<-mydata[!ERCC_ind,N1_1]
N1_2_DAT<-mydata[!ERCC_ind,N1_2]
N1_3_DAT<-mydata[!ERCC_ind,N1_3]

N2_1_DAT<-mydata[!ERCC_ind,N2_1]
N2_2_DAT<-mydata[!ERCC_ind,N2_2]
N2_3_DAT<-mydata[!ERCC_ind,N2_3]

N3_1_DAT<-mydata[!ERCC_ind,N3_1]
N3_2_DAT<-mydata[!ERCC_ind,N3_2]
N3_3_DAT<-mydata[!ERCC_ind,N3_3]


N123_123_DAT<-cbind(N1_1_DAT,N1_3_DAT,N2_1_DAT,N2_2_DAT,N2_3_DAT,N3_1_DAT,N3_2_DAT,N3_3_DAT)
dim(N123_123_DAT)


#Individual NA19098
N1_DAT<-cbind(N1_1_DAT,N1_3_DAT)
#Individual NA19101
N2_DAT<-cbind(N2_1_DAT,N2_2_DAT,N2_3_DAT)
#Individual NA19239
N3_DAT<-cbind(N3_1_DAT,N3_2_DAT,N3_3_DAT)

#Individuals' labels
LABEL_INDIVIDUAL<-c(rep("NA19098",dim(N1_DAT)[2]),rep("NA19101",dim(N2_DAT)[2]),rep("NA19239",dim(N3_DAT)[2]))

#Batch labels
LABEL_REP<-c(rep(1,dim(N1_1_DAT)[2]),rep(3,dim(N1_3_DAT)[2]),rep(1,dim(N2_1_DAT)[2]),rep(2,dim(N2_2_DAT)[2]),rep(3,dim(N2_3_DAT)[2]),rep(1,dim(N3_1_DAT)[2]),rep(2,dim(N3_2_DAT)[2]),rep(3,dim(N3_3_DAT)[2]))






######Using ERCC spike-ins to estimate BETA in Tung study#########

#### Code was adapted from https://github.com/jdblischak/singleCellSeq/blob/master/analysis/capture-efficiency.Rmd######


anno <- read.table("E:/RNAseqProject/tung2017batch/annotation.txt", header = TRUE,stringsAsFactors = FALSE)
head(anno)


ercc <- read.table("E:/RNAseqProject/tung2017batch/ercc-info.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
colnames(ercc) <- c("num", "id", "subgroup", "conc_mix1", "conc_mix2", "expected_fc", "log2_mix1_mix2")
head(ercc)
stopifnot(nrow(ercc) == 92)

molecules <- read.table("E:/RNAseqProject/tung2017batch/molecules.txt", header = TRUE,stringsAsFactors = FALSE)
quality_single_cells <- scan("E:/RNAseqProject/tung2017batch/quality-single-cells.txt",what = "character")

qc <- read.table("E:/RNAseqProject/tung2017batch/qc-ipsc.txt", header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(qc) == sum(anno$well != "bulk"))
length(molecules$NA19098.r1.A03[rownames(N1_DAT)])



#Prepare single cell molecule data
molecules_single <- molecules[, anno$well != "bulk"]
anno_single <- anno[anno$well != "bulk", ]
#Remove genes with zero read counts in the single cells.
expressed_single <- rowSums(molecules_single) > 0
molecules_single <- molecules_single[expressed_single, ]
dim(molecules_single)

overexpressed_genes <- rownames(molecules_single)[apply(molecules_single, 1,function(x) any(x >= 1024))]
length(overexpressed_genes)

#####Calculate number of ERCC molecules added to each well
summary(ercc$conc_mix1)

# Dilute 1:2500
ercc_conc_diluted <- ercc$conc_mix1 / 2500
# Dilute 1:20
ercc_conc_lysis <- ercc_conc_diluted / 20
ercc_molecules_lysis <- ercc_conc_lysis *
    20 * # Number of uL of lysis buffer
    1/10^18 * # Number of attomoles in a mole
    6.02214179e23 # Number of molecules in a mole
# 9 uL added to chip
ercc_molecules_chip <- ercc_molecules_lysis * 9 / 20
summary(ercc_molecules_chip)
# 9 nL per well
ercc_molecules_well <- ercc_molecules_lysis * 9e-3 / 20
summary(ercc_molecules_well)
sum(ercc_molecules_well)
sum(ercc_molecules_well >= 1)
sum(ercc_molecules_well > 1024)
sum(ercc_molecules_well %% 2 == 0)

ercc_index <- grep("ERCC", rownames(molecules_single))
length(ercc_index)
#####capture efficiency#########
expected_ercc_molecules <- read.table("E:/RNAseqProject/tung2017batch/expected-ercc-molecules.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)

efficiency <- numeric(length = ncol(molecules_single))
total_ercc_molecules <- sum(ercc_molecules_well)
for (i in 1:ncol(molecules_single)) {
    efficiency[i] <- sum(molecules_single[ercc_index, i]) / total_ercc_molecules
}
summary(efficiency)
names(efficiency)<-colnames(molecules_single)



####Store raw data for Tung for further study#######
save.image(file="E:/RNAseqProject/Illustrator_bayNorm/bayNorm_papercode/RealData/Tung_study/Load_Tung.RData")
