pathtttt<-"E:/RNAseqProject/tung2017batch/MMEADJ_ERCCBETA/array_batchcheck_saversample_tr_default"

load(paste(pathtttt,"/MT_N2_12.RData",sep=''))
load(paste(pathtttt,"/MT_N2_13.RData",sep=''))
load(paste(pathtttt,"/MT_N2_23.RData",sep=''))

load(paste(pathtttt,"/MT_N3_12.RData",sep=''))
load(paste(pathtttt,"/MT_N3_13.RData",sep=''))
load(paste(pathtttt,"/MT_N3_23.RData",sep=''))

load(paste(pathtttt,"/MT_N1_13.RData",sep=''))



#save EPS files into other file
pathhhh<-"E:/RNAseqProject/tung2017batch/MMEADJ_ERCCBETA/mean_bayNorm/meanBATCH_CHECK/forpaper/array_saversample_default"
path_fun<-function(filename){
    qq<-paste(pathhhh,filename,sep='')
    return(qq)
}

source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")
methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC','DCA')






######begin multiple plots######
textsize=10
geom_text_size=1.5
legend_key_size=1
#replot
BAR_M_N1_13<-BATCH_BAR_FUN2(M_list_N1_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)


BAR_M_N2_13<-BATCH_BAR_FUN2(M_list_N2_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N2_12<-BATCH_BAR_FUN2(M_list_N2_12,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N2_23<-BATCH_BAR_FUN2(M_list_N2_23,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)


BAR_M_N3_13<-BATCH_BAR_FUN2(M_list_N3_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N3_12<-BATCH_BAR_FUN2(M_list_N3_12,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N3_23<-BATCH_BAR_FUN2(M_list_N3_23,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)




library(gridExtra)
library(ggpubr)
library(cowplot)

qq<-plot_grid(BAR_M_N2_13 + theme(legend.position="none"),BAR_M_N2_12 + theme(legend.position="none"),BAR_M_N2_23 + theme(legend.position="none"), BAR_M_N3_13 + theme(legend.position="none"),BAR_M_N3_12 + theme(legend.position="none"),BAR_M_N3_23 + theme(legend.position="none"),BAR_M_N1_13 + theme(legend.position="none"),ncol=3,nrow=3)+ draw_grob(get_legend(BAR_M_N2_13), (1.5)/3, 0, 1/3, 0.5)
qq

ggsave(filename=path_fun('/Tung_Mbatch_ERCCBETA_array_tr_default.pdf'),plot=qq, width =8.2, height =9.5,units='in',device='pdf')


######################individual plots##########
source("E:/RNAseqProject/tung2017batch/FINAL/BATCH_DE_CHECK/BATCH_BAR_FUN.R")
methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC')
textsize=10
geom_text_size=1.8
legend_key_size=0.3
#replot
BAR_M_N1_13<-BATCH_BAR_FUN2(M_list_N1_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N2_13<-BATCH_BAR_FUN2(M_list_N2_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N2_12<-BATCH_BAR_FUN2(M_list_N2_12,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N2_23<-BATCH_BAR_FUN2(M_list_N2_23,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)

BAR_M_N3_13<-BATCH_BAR_FUN2(M_list_N3_13,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)
BAR_M_N3_12<-BATCH_BAR_FUN2(M_list_N3_12,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)
BAR_M_N3_23<-BATCH_BAR_FUN2(M_list_N3_23,methodnames,textsize=textsize,geom_text_size=geom_text_size,legend_key_size=legend_key_size)


#width =2.75, height =3.8
#MAST###### 
width<-4.5
height=3.8

###M_N1
#postscript(path_fun('/ERCC_BETA/M_N1_13.eps'), width =width, height =height)
pdf(path_fun('/ERCC_BETA/M_N1_13.pdf'), width =width, height =height)
BAR_M_N1_13
dev.off()

###M_N2
postscript(path_fun('/ERCC_BETA/M_N2_13.eps'), width =width, height =height)
BAR_M_N2_13
dev.off()

postscript(path_fun('/ERCC_BETA/M_N2_12.eps'), width =width, height =height)
BAR_M_N2_12
dev.off()

postscript(path_fun('/ERCC_BETA/M_N2_23.eps'), width =width, height =height)
BAR_M_N2_23
dev.off()

###M_N3
postscript(path_fun('/ERCC_BETA/M_N3_13.eps'), width =width, height =height)
BAR_M_N3_13
dev.off()

postscript(path_fun('/ERCC_BETA/M_N3_12.eps'), width =width, height =height)
BAR_M_N3_12
dev.off()

postscript(path_fun('/ERCC_BETA/M_N3_23.eps'), width =width, height =height)
BAR_M_N3_23
dev.off()



###M_AVG####
#####mean of false positive rates across pair of batches:######
BAR_M_N1_13$data$`Number of detected DE genes`


BAR_M_N1_13$data$`Normalization method`
qq<-as.data.frame(BAR_M_N1_13$data)


xxx_data<-cbind(BAR_M_N1_13$data,BAR_M_N2_13$data,BAR_M_N2_12$data,BAR_M_N2_23$data,BAR_M_N3_13$data,BAR_M_N3_12$data,BAR_M_N3_23$data)
dim(xxx_data)

num_gene<-length(M_list_N1_13[[1]])
kkk_data<-data.frame(apply(xxx_data[,c(1,4,7,10,13,16,19)],1,mean)/num_gene,qq$`Adjusted P-values threshold`,qq$`Normalization method`)

geom_text_size<-1.2
textsize<-8
geom_text_size=1.8
legend_key_size=0.3


methodnames<-c('bayNorm','SCnorm','SAVER','scImpute','Scaling','MAGIC','DCA')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette )<-c('NULL','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','DCA')
cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% methodnames)]

colnames(kkk_data)<-c('Avg_FP','AdjPval','Methods')
kkk_data_bar<-ggplot(data=kkk_data, aes(x=kkk_data[,2], y=kkk_data[,1], fill=kkk_data[,3])) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    #geom_text(aes(label=round(kkk_data[,1],4)), vjust=1.6, color="black", position = position_dodge(0.9), size=geom_text_size)+
    labs(x = 'Adjusted P values threshold',y='Averaged False positive rates',fill='Normalization methods')+
    ggtitle("") +
    #scale_fill_brewer(palette="Paired")+
    scale_fill_manual(values=cbbPalette2)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size =textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.position="top",legend.key.size = unit(legend_key_size,"line"))

avg_005<-kkk_data$Avg_FP[kkk_data$AdjPval==0.05]
names(avg_005)<-kkk_data$Methods[kkk_data$AdjPval==0.05]
avg_001<-kkk_data$Avg_FP[kkk_data$AdjPval==0.01]
names(avg_001)<-kkk_data$Methods[kkk_data$AdjPval==0.01]
avg_01<-kkk_data$Avg_FP[kkk_data$AdjPval==0.1]
names(avg_01)<-kkk_data$Methods[kkk_data$AdjPval==0.1]

save(avg_005,avg_001,avg_01,file="E:/RNAseqProject/tung2017batch/MMEADJ_ERCCBETA/avg_vectors_default.RData")

kkk_data_bar

width<-4.5
height=3.8
g <- ggplot_build(kkk_data_bar)
unique(g$data[[1]]["fill"])

# ccc<-rev(c("#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3"))
# plot(seq(1,6),col=rev(ccc),pch=16,cex=2)



pdf(path_fun('/M_AVG_array_tr_BETA.pdf'), width =width, height =height)
kkk_data_bar
dev.off()
