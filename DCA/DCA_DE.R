source("E:/RNAseqProject/MANY_DE_FUN.R")

##dca ISLAM
load("E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_Islam/DCA_Islam.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_Islam,NumCells = c(48,44))
save(DCA_Islam,M_DCA,file="E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_Islam/DCA_Islam.RData")


##dca H1
load("E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_H1/DCA_H1.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_H1,NumCells = c(92,92))
save(DCA_H1,M_DCA,file="E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_H1/DCA_H1.RData")
load("E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_H1/DCA_H1.RData")

##dca H9
load("E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_H9/DCA_H9.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_H9,NumCells = c(91,91))
save(DCA_H9,M_DCA,file="E:/RNAseqProject/RAW_REAL/DCA_norm/DCA_H9/DCA_H9.RData")



#SIM 1
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_1.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_1.RData")

#SIM 01_005
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_01_005.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_01_005.RData")

#SIM 005_01
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_005_01.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_005_01.RData")


#SIM 005_005
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_005_005.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_005_005.RData")


#SIM noDE_01_005
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_noDE_01_005.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_noDE_01_005.RData")

length(which(M_DCA$adjpval<0.05))
plot(rowMeans(DCA_out[,seq(1,100)]),rowMeans(DCA_out[,-seq(1,100)]),log='xy')
abline(0,1)

#SIM noDE_01_005_H1
load("E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_noDE_01_005_H1.RData")
M_DCA<-SCnorm_runMAST3(Data=DCA_out,NumCells = c(100,100))
save(DCA_out,M_DCA,file="E:/RNAseqProject/RAW_6SIM/DCA_SIM/DCA_SIM_noDE_01_005_H1.RData")








