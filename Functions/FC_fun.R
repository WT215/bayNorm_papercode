

FC_fun<-function(Inputdat,CONDITION,DATA_list,textsize=1,legend.key.size=1,colourval=NULL){

MedExp <- log(apply(Inputdat, 1, function(x) median(x[x != 0])))
# split into 4 equally sized groups:
grpnum <- 6
splitby <- sort(MedExp)
grps <- length(splitby)/grpnum
sreg <- split(splitby, ceiling(seq_along(splitby)/grps))


#######begin many plots########



#######Figure 2a############

#logfold='My'
logfold='SCnorm'
 if (logfold=='SCnorm')
{
    
    library(foreach)
     
     fc_list<-foreach(i=1:length(DATA_list))%do%{
         

             fc <- apply(DATA_list[[i]][,CONDITION==1], 1, function(x) mean((x[x!=0]))) / apply(DATA_list[[i]][,CONDITION==2], 1, function(x) mean((x[x!=0])))
         #fc <- apply(DATA_list[[i]][,CONDITION==1], 1, function(x) mean((x))) / apply(DATA_list[[i]][,CONDITION==2], 1, function(x) mean((x)))
             fc<-log2(fc)

         return(fc)
         
     }
     
   names(fc_list)<-names(DATA_list)
}


#fc_list<-fc_list[-2]

library(foreach)
Bar_data<-foreach(i=1:length(sreg),.combine=rbind)%:%
    
    foreach(j = 1:length(fc_list),.combine=rbind)%do%{
        qq<-cbind( fc_list[[j]][names(sreg[[i]])],rep(i,length(fc_list[[j]][names(sreg[[i]])])),rep(names(fc_list)[j],length(fc_list[[j]][names(sreg[[i]])])))
        
        return(qq)
    }
Bar_data<-as.data.frame(Bar_data)

colnames(Bar_data)<-c('fc','Group','Method')
Bar_data$fc<-as.numeric(as.character(Bar_data$fc))

Bar_data$Group<-factor(Bar_data$Group,levels=unique(Bar_data$Group))
Bar_data$Method<-factor(Bar_data$Method,levels=unique(Bar_data$Method))




library(ggplot2)

FC_H1<-ggplot(data=Bar_data, aes(x=Bar_data$Group, y=Bar_data$fc, fill=Bar_data$Method)) +
    geom_boxplot(size=0.05,outlier.size = 0.001)+
    labs(x = "Expression group",y='Log2 Fold-Change',fill='Normalization methods')+ggtitle("") +
    geom_hline(aes(yintercept=0),size=0.15)+
    #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=colourval)+
  
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(legend.key.size,"line"),legend.position ='top')



return(FC_H1)
}
