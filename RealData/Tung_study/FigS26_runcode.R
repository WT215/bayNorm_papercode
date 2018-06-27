source("E:/RNAseqProject/MANY_DE_FUN.R")
source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")

load(file='Tung_norms.RData')
pathtttt<-"E:/RNAseqProject/tung2017batch/MMEADJ_ERCCBETA/array_batchcheck_saversample_tr_default"


names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N2<-LABEL_REP[colnames(N2_DAT)]

N2_meanbay<-Bay_2
N2_scnorm<-scnorm_norm[,colnames(N2_DAT)]
N2_1oversaver<-saver_array_N2
N2_scImpute<-scImpute_N2_DAT
N2_RB<-RB_N2
N2_MAGIC<-MAGIC_TUNG[,colnames(N2_DAT)]
N2_DCA<-DCA_tungall[,colnames(N2_DAT)]

methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC','DCA')

N2_LIST<-list(N2_meanbay,N2_scnorm,N2_1oversaver,N2_scImpute,N2_RB,N2_MAGIC,N2_DCA)


names(N2_LIST)<-methodnames



#N2(individual NA19101): 1 2 #######
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N2<-LABEL_REP[colnames(N2_DAT)]
q_temp<-BATCH_N2[BATCH_N2!=3]



DATA_list_N2_12<-foreach(i=1:length(N2_LIST))%do%{
    if(length(dim(N2_LIST[[i]]))==3){
        return(N2_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N2_LIST[[i]]))==2){
        return(N2_LIST[[i]][,names(q_temp)])
    }
}


M_list_N2_12<-M_TEST_FUN(DATA_list_N2_12,grs=table(q_temp),methodnames=methodnames)


BAR_M_N2_12<-BATCH_BAR_FUN(M_list_N2_12,methodnames=methodnames)
#BAR_T_N2_12<-BATCH_BAR_FUN(T_list_N2_12,methodnames=methodnames)



save(M_list_N2_12,BAR_M_N2_12,file=paste(pathtttt,"/MT_N2_12.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N2_12.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N2_12
dev.off()

#N2(individual NA19101): 1 3#####
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N2<-LABEL_REP[colnames(N2_DAT)]
q_temp<-BATCH_N2[BATCH_N2!=2]


DATA_list_N2_13<-foreach(i=1:length(N2_LIST))%do%{
    if(length(dim(N2_LIST[[i]]))==3){
        return(N2_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N2_LIST[[i]]))==2){
        return(N2_LIST[[i]][,names(q_temp)])
    }
}


M_list_N2_13<-M_TEST_FUN(DATA_list_N2_13,grs=table(q_temp),methodnames)

BAR_M_N2_13<-BATCH_BAR_FUN(M_list_N2_13,methodnames)


save(M_list_N2_13,BAR_M_N2_13,file=paste(pathtttt,"/MT_N2_13.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N2_13.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N2_13
dev.off()


#N2(individual NA19101): 2 3########
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N2<-LABEL_REP[colnames(N2_DAT)]
q_temp<-BATCH_N2[BATCH_N2!=1]


DATA_list_N2_23<-foreach(i=1:length(N2_LIST))%do%{
    if(length(dim(N2_LIST[[i]]))==3){
        return(N2_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N2_LIST[[i]]))==2){
        return(N2_LIST[[i]][,names(q_temp)])
    }
}

M_list_N2_23<-M_TEST_FUN(DATA_list_N2_23,grs=table(q_temp),methodnames)


BAR_M_N2_23<-BATCH_BAR_FUN(M_list_N2_23,methodnames)


save(M_list_N2_23,BAR_M_N2_23,file=paste(pathtttt,"/MT_N2_23.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N2_23.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N2_23
dev.off()


######N3(individual NA19239)########
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N3<-LABEL_REP[colnames(N3_DAT)]

N3_meanbay<-Bay_3
N3_scnorm<-scnorm_norm[,colnames(N3_DAT)]
#N3_1oversaver<-saver_N3$estimate
N3_1oversaver<-saver_array_N3
N3_scImpute<-scImpute_N3_DAT
N3_RB<-RB_N3
N3_MAGIC<-MAGIC_TUNG[,colnames(N3_DAT)]
N3_DCA<-DCA_tungall[,colnames(N3_DAT)]


methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC','DCA')

N3_LIST<-list(N3_meanbay,N3_scnorm,N3_1oversaver,N3_scImpute,N3_RB,N3_MAGIC,N3_DCA)
names(N3_LIST)<-methodnames




#N3(individual NA19239): 1 2#######
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N3<-LABEL_REP[colnames(N3_DAT)]
q_temp<-BATCH_N3[BATCH_N3!=3]


DATA_list_N3_12<-foreach(i=1:length(N3_LIST))%do%{
    if(length(dim(N3_LIST[[i]]))==3){
        return(N3_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N3_LIST[[i]]))==2){
        return(N3_LIST[[i]][,names(q_temp)])
    }
}


M_list_N3_12<-M_TEST_FUN(DATA_list_N3_12,grs=table(q_temp),methodnames)



BAR_M_N3_12<-BATCH_BAR_FUN(M_list_N3_12,methodnames)

save(M_list_N3_12,BAR_M_N3_12,file=paste(pathtttt,"/MT_N3_12.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N3_12.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N3_12
dev.off()

#N3(individual NA19239): 1 3#####
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N3<-LABEL_REP[colnames(N3_DAT)]
q_temp<-BATCH_N3[BATCH_N3!=2]


DATA_list_N3_13<-foreach(i=1:length(N3_LIST))%do%{
    if(length(dim(N3_LIST[[i]]))==3){
        return(N3_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N3_LIST[[i]]))==2){
        return(N3_LIST[[i]][,names(q_temp)])
    }
}


M_list_N3_13<-M_TEST_FUN(DATA_list_N3_13,grs=table(q_temp),methodnames)


BAR_M_N3_13<-BATCH_BAR_FUN(M_list_N3_13,methodnames)

save(M_list_N3_13,BAR_M_N3_13,file=paste(pathtttt,"/MT_N3_13.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N3_13.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N3_13
dev.off()



#N3(individual NA19239): 2 3########
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N2<-LABEL_REP[colnames(N3_DAT)]
q_temp<-BATCH_N3[BATCH_N3!=1]


DATA_list_N3_23<-foreach(i=1:length(N3_LIST))%do%{
    if(length(dim(N3_LIST[[i]]))==3){
        return(N3_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N3_LIST[[i]]))==2){
        return(N3_LIST[[i]][,names(q_temp)])
    }
}

M_list_N3_23<-M_TEST_FUN(DATA_list_N3_23,grs=table(q_temp),methodnames)



BAR_M_N3_23<-BATCH_BAR_FUN(M_list_N3_23,methodnames)


save(M_list_N3_23,BAR_M_N3_23,file=paste(pathtttt,"/MT_N3_23.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N3_23.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N3_23
dev.off()




#####N1 (individual NA19098)########
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N1<-LABEL_REP[colnames(N1_DAT)]

N1_meanbay<-Bay_1
N1_scnorm<-scnorm_norm[,colnames(N1_DAT)]
#N1_1oversaver<-saver_N1$estimate
N1_1oversaver<-saver_array_N1
N1_scImpute<-scImpute_N1_DAT
N1_RB<-RB_N1
N1_MAGIC<-MAGIC_TUNG[,colnames(N1_DAT)]
N1_DCA<-DCA_tungall[,colnames(N1_DAT)]


methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC','DCA')

N1_LIST<-list(N1_meanbay,N1_scnorm,N1_1oversaver,N1_scImpute,N1_RB,N1_MAGIC,N1_DCA)
names(N1_LIST)<-methodnames


source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")

#N1 (individual NA19098): 1 3#######
names(LABEL_REP)<-colnames(cbind(N1_DAT,N2_DAT,N3_DAT))
BATCH_N1<-LABEL_REP[colnames(N1_DAT)]
q_temp<-BATCH_N1[BATCH_N1!=2]


DATA_list_N1_13<-foreach(i=1:length(N1_LIST))%do%{
    if(length(dim(N1_LIST[[i]]))==3){
        return(N1_LIST[[i]][,names(q_temp),])
    } else if(length(dim(N1_LIST[[i]]))==2){
        return(N1_LIST[[i]][,names(q_temp)])
    }
}

M_list_N1_13<-M_TEST_FUN(DATA_list_N1_13,grs=table(q_temp),methodnames=methodnames)

BAR_M_N1_13<-BATCH_BAR_FUN(M_list_N1_13,methodnames=methodnames)


save(M_list_N1_13,BAR_M_N1_13,file=paste(pathtttt,"/MT_N1_13.RData",sep=''))
jpeg(paste(pathtttt,"/BAR_M_N1_13.jpeg",sep=''), width = 40, height = 30, units = 'in', res = 300)
BAR_M_N1_13
dev.off()










