#Benchmark DE genes and raw data were kindly provided by Jaakkola (doi: 10.1093/bib/bbw057)
Top1000_Jaakkola<-read.table(file="~/Top1000_Jaakkola.txt")
DAT_Jaakkola<-read.delim("~/Islam_overlap.txt", sep = "\t" , header = T)

length(intersect(Top1000_Jaakkola$V1,rownames(DAT_Jaakkola)))
DAT_Jaakkola<-as.matrix(DAT_Jaakkola)
length(which(rowSums(DAT_Jaakkola)==0))
drop<-which(rowSums(DAT_Jaakkola)==0)

TRUE_LABEL<-rep(0,dim(DAT_Jaakkola)[1])
names(TRUE_LABEL)<-rownames(DAT_Jaakkola)
TRUE_LABEL[which(rownames(DAT_Jaakkola) %in% as.character(Top1000_Jaakkola$V1))]<-1
grs = table(substr(colnames(DAT_Jaakkola),1,2))
DAT_Jaakkola<-DAT_Jaakkola[-drop,]

#TRUE_LABEL is the label of DE genes, used as benchmark
TRUE_LABEL<-TRUE_LABEL[-drop]









save('Load_Islam.RData')



