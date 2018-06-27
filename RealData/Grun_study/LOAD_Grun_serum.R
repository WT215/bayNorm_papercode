#####import raw data########
library(data.table)

data <- fread("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/GSE54695_data_transcript_counts.txt")

Gene.Ids = data$GENENAME
Cell.Ids = names(data)[-1]

RawCounts = cbind(subset(data, select = Cell.Ids[grep("SC_serum", Cell.Ids)]),subset(data, select = Cell.Ids[grep("RNA_serum", Cell.Ids)]))

Cell.Colour = c( rep("lightpink3",length(grep("SC_serum", Cell.Ids))),rep("darkolivegreen3",length(grep("RNA_serum", Cell.Ids))) )

Batch = c(rep(1, 40), rep(2, 40), rep(1, 40), rep(2, 40))
dim(RawCounts)

#Transforming the data into UMI counts
# Function provided by Jong Kyoung Kim (EMBL-EBI)
UMICount <- function(MoleculeCount, UMILength)
{
  # MoleculeCount is the normalized count
  M = 4^UMILength
  UMICount = M*(1-exp(-MoleculeCount/M))
  return(UMICount)
}
CountsUMI = round(UMICount(RawCounts, 4))


Pou5f1.per.cell <- as.numeric(CountsUMI[which(Gene.Ids == "Pou5f1"),])
counts.per.cell <- colSums(CountsUMI)
genes.per.cell <- apply(CountsUMI, 2, function(x) sum( x>0 ))
ercc.per.cell <- colSums(CountsUMI[grep("ERCC", Gene.Ids),])

###After filtering
CountsUMI_1 <- CountsUMI[,Pou5f1.per.cell >= 10, with = FALSE]
Cell.Colour_1 <- Cell.Colour[Pou5f1.per.cell >= 10]
Batch_1 <- Batch[Pou5f1.per.cell >= 10]
table(Cell.Colour_1)
dim(CountsUMI_1)
CountsUMI_2 = CountsUMI_1[rowSums(CountsUMI_1) >= 50, ]
Gene.Ids_2 = Gene.Ids[rowSums(CountsUMI_1) >= 50]
dim(CountsUMI_2)

#Spike-in genes information#######
SpikeInfo<-fread("E:/RNAseqProject/NEWPROJECT_PAPERS/Validation of noise models for single-cell transcriptomics/cms_095046.txt")
names(SpikeInfo)[2]<-"ERCC_ID"
names(SpikeInfo)[4]<-"concentration_in_Mix_1"


SpikeInfo <- SpikeInfo[SpikeInfo$ERCC_ID %in% Gene.Ids_2[grep("ERCC", Gene.Ids_2)],]
SpikeInfo$MoleculesPerCell <- SpikeInfo$concentration_in_Mix_1 * (1e-18) *(6.022e23) * (1/2500000)
# To confirm with the 3.3% capture indicated by Grun et al
SpikeOut <- data.table("ERCC_ID" = Gene.Ids_2[grep("ERCC", Gene.Ids_2)],"ERCC_MeanCount" = rowMeans(CountsUMI_2[grep("ERCC", Gene.Ids_2),]))
SpikeOut = merge(SpikeInfo, SpikeOut, by = "ERCC_ID")
SpikeInfoFilter = subset(SpikeInfo, select = c(ERCC_ID, MoleculesPerCell))
Tech = grepl("ERCC", Gene.Ids_2)
Total_ERCC<-sum(SpikeInfoFilter$MoleculesPerCell)
BETA_ERCC<-colSums(CountsUMI_2[Tech,])/Total_ERCC
summary(SpikeOut$ERCC_MeanCount/SpikeOut$MoleculesPerCell)


#Separating expression counts for each condition
CountsUMI_SC = as.matrix(CountsUMI_2)[, grep("SC_serum",colnames(CountsUMI_2))]
CountsUMI_RNA = as.matrix(CountsUMI_2)[, grep("RNA_serum",colnames(CountsUMI_2))]
rownames(CountsUMI_SC) <- Gene.Ids_2
rownames(CountsUMI_RNA) <- Gene.Ids_2

CountsUMI_SC_ERCC<-CountsUMI_SC[Tech,]
CountsUMI_RNA_ERCC<-CountsUMI_RNA[Tech,]
CountsUMI_SC<-CountsUMI_SC[-Tech,]
CountsUMI_RNA<-CountsUMI_RNA[-Tech,]

CountsUMI_SC_serum<-CountsUMI_SC


save.image('Grun_2014_RAW_serum.RData')