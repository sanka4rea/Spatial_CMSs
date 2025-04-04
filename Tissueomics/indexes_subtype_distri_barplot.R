# return plot and wilcoxn p values
distri_barplot <- function(Matrix = EQ_index_fm, index = EQ_index_fm$lym, wilcox_pair=NULL,title=NULL,labels=EQ_index_fm$Subtype,
                           colors=c("#fc6c85","#c364c5","#ffa343","#1cac78","#1974d2")){ # c("#E69E00","#0070B0","#CA78A6","#009C73")
  # "#00468B","#ED0000","#42B540","#0099B4","dimgray"
  image_str<- data.frame(index=index) 
  image_str$cms <- labels
  image_str$color <- labels
  if(is.null(wilcox_pair)){
    P_wilcox <-NA
  }else{
    P_wilcox <- wilcox.test(image_str$index[which(image_str$cms==as.character(wilcox_pair))],
                            image_str$index[which(image_str$cms!=as.character(wilcox_pair))])
  }
  strtum_mean <- lapply(sort(unique(as.numeric(labels))), function(xx){
    tmp <- mean(image_str$index[which(image_str$cms==as.character(xx))],na.rm=TRUE)
    return(tmp)
  })
  strtum_mean <- as.vector(do.call(cbind,strtum_mean))
  
  strtum_std <- lapply(sort(unique(as.numeric(labels))), function(xx){
    tmp <- sd(image_str$index[which(image_str$cms==as.character(xx))],na.rm=TRUE)/sqrt(length(which(image_str$cms==as.character(xx))))
    return(tmp)
  })
  strtum_std <- as.vector(do.call(cbind,strtum_std))

  strtum <- cbind(strtum_mean,strtum_std)
  strtum <- as.data.frame(strtum)
  strtum$CMS <- as.character(sort(unique(as.numeric(labels))) ) 
  colnames(strtum) <- c("Average","SD","CMS")
  lymPlot <-ggplot(strtum,aes(x=CMS, y=Average,fill=CMS,color=CMS)) + 
    geom_bar(stat = "identity",position = "dodge",colour="black",cex=1,width = 0.6)+ #
    #scale_fill_brewer(palette = "Pastel1")+
    labs(title=title, x="", y="Ratio (%)")+
    theme(axis.title.x = element_text(size = 20,family="sans"))+
    theme(axis.title.y = element_text(size = 20,family="sans"))+ #, angle = 90
    scale_fill_manual(values=colors)+ 
    scale_color_manual(values=colors)+
    theme(axis.text.x = element_text(size=20,family="sans"),axis.text.y = element_text(size=20,family="sans"),
          plot.title = element_text( color="black", size=20,face="plain",hjust = 0.5,family ="sans"))+
    # theme(aspect.ratio = 1.3) +
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
    )+ # axis.ticks = element_blank()
    geom_errorbar(aes(ymin=Average-SD, ymax=Average+SD), width=.2,color="black",
                  position=position_dodge(.8),cex=1)
  # +coord_cartesian(ylim=c(0.01, 0.042))
  return(list(lymPlot,P_wilcox))
}



