####sup fun###

qwerfun<-function(DenDat,linesize=1,CAPTION=''){
  
  pp<-unique(DenDat$`Normalization method`)
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  names(cbbPalette )<-c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','Raw')
  cbbPalette2<-cbbPalette[which(names(cbbPalette) %in% pp)]
  
  
  Gg<-unique(DenDat$Genename)
  typp<-unique(DenDat$Condition)
  
  filenn<-paste('Gene: ',Gg,sep='')
  
  density_plot<-ggplot() +
    #<-ggplot(DenDat, aes(x=DenDat$`Normalized count`)) +
    #geom_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,position="identity",show.legend=FALSE)+
    stat_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,geom="line",position="identity")+
    
    
    
    scale_colour_manual( values = cbbPalette2)+
    #add xlim for other than Grun data
    xlim(0, max(DenDat$`Normalized count`[which(DenDat$`Normalization method`=='smFISH')]))+
    labs(x = "Gene expression",y='density',color='Datasets')+
    ggtitle(filenn,subtitle=CAPTION)+
    #guides (colour = guide_legend (override.aes=list(colour=NA)))+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1,"line"),legend.position = c(0.8,0.9),legend.background=element_blank(),legend.key =element_blank())
  
  
  
  g1 <- ggplotGrob(density_plot)
  # Check out the grobs
  library(grid)
  #grid.ls(grid.force(g1))
  # Get names of 'label' grobs.
  names.grobs <- grid.ls(grid.force(g1), print=F)$name 
  labels <- names.grobs[which(grepl("^label", names.grobs))]
  
  for(i in seq_along(labels)) {
    g1 <- editGrob(grid.force(g1), gPath(labels[i]), grep = TRUE,  gp = gpar(col =cbbPalette2[i]))
  }
  
  
  density_plot_v2<-as_ggplot(g1)
  
  #return(density_plot)
  return(density_plot)
}



qwerfun_ind<-function(DenDat,filenn='',linesize=1,CAPTION=''){
  
  pp<-unique(DenDat$`Normalization method`)
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  names(cbbPalette )<-c('smFISH','bayNorm','SCnorm','Scaling','SAVER','scImpute','MAGIC','Raw')
  cbbPalette2<-cbbPalette[match(pp,names(cbbPalette))]

  
  Gg<-unique(DenDat$Genename)
  typp<-unique(DenDat$Condition)
  
  filenn<-paste('Gene: ',Gg,sep='')
  
  density_plot<-ggplot() +
    #<-ggplot(DenDat, aes(x=DenDat$`Normalized count`)) +
    #geom_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),bw=bw,size=linesize,position="identity",show.legend=FALSE)+
    stat_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,geom="line",position="identity")+
    scale_colour_manual( values = cbbPalette2)+
    #labs(x = "Normalized count",y='density',color='Normalization methods')+
    labs(x = "Normalized count",y='density',color='')+
    ggtitle(filenn,subtitle=CAPTION)+
    guides (colour = guide_legend (override.aes=list(colour=NA)))+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1.5,"line"),legend.position = c(0.75,0.85),legend.background=element_blank(),legend.key =element_blank())
  
  g1 <- ggplotGrob(density_plot)
  # Check out the grobs
  library(grid)
  #grid.ls(grid.force(g1))
  # Get names of 'label' grobs.
  names.grobs <- grid.ls(grid.force(g1), print=F)$name 
  labels <- names.grobs[which(grepl("^label", names.grobs))]
  
  for(i in seq_along(labels)) {
    g1 <- editGrob(grid.force(g1), gPath(labels[i]), grep = TRUE,  gp = gpar(col =cbbPalette2[i]))
  }
  
  # Draw it
  #grid.newpage()
  #density_plot_v2<- grid.draw(g1, recording=F)
  density_plot_v2<-as_ggplot(g1)
  
  return( density_plot_v2)
}


Density_fun<-function(list_vec,methodd,Gg,filenn='',textsize=1,linesize=1,CAPTION='',bw=10){
  
  library(ggplot2)
  library(foreach)
  DenDat<-foreach(i=1:length(list_vec),.combine=rbind)%do%{
    aa<-list_vec[[i]]
    bb<-cbind(aa,rep(methodd[i],length(aa)),rep(Gg,length(aa)))
    return(bb)
  }
  colnames(DenDat)<-c('Normalized count','Normalization method','Genename')

  DenDat<-as.data.frame(DenDat)
  DenDat$`Normalized count`<-as.numeric(as.character(DenDat$`Normalized count`))
  DenDat$`Normalization method`<-factor(DenDat$`Normalization method`,levels=unique(DenDat$`Normalization method`))
  
  cbPalette <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  filenn<-paste('Gene: ',Gg,sep='')
  
  density_plot<-ggplot() +
    stat_density(aes(x=DenDat$`Normalized count`,group=DenDat$`Normalization method`,colour=DenDat$`Normalization method`),size=linesize,bw=bw,geom="line",position="identity")+
    scale_colour_manual( values = cbPalette)+
    labs(x = "Normalized count",y='density',color='Normalization methods',caption=CAPTION)+
    ggtitle(filenn)+
    theme(legend.text = element_text(size = textsize),legend.title  = element_text(size = textsize),plot.title = element_text(size = textsize),axis.title = element_text(size = textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.subtitle = element_text(size = textsize),plot.caption =  element_text(size = textsize),axis.text=element_text(size=textsize) ,legend.key.size = unit(1.5,"line"),legend.position = c(0.85,0.85))
  return(list(density_plot=density_plot,DenDat=DenDat))
}

