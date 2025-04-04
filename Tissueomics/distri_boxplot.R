# return plot and wilcoxn p values
distri_barplot <- function(Matrix = Im_erqiClinical_1$Entropy, labels = Im_erqiClinical_1$snf4more, title=NULL,ylab=NULL,
                           colors=c("#00468B","#ED0000","#42B540","#0099B4","dimgray")){ # c("#E69E00","#0070B0","#CA78A6","#009C73")
  
  dck13 <- data.frame(dis=Matrix,label=labels)
  plot <- ggplot(dck13, aes(x=label, y=dis,fill=label)) +  #,color=label
    geom_boxplot(notch=F,notchwidth = 0.1,width = 0.4,varwidth = FALSE)+ #lwd=1.5,
    stat_boxplot(geom ='errorbar',width = 0.3)+
    labs(title=title, x="", y=ylab)+
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text( color="black", size=15,face="plain",hjust = 0.5),
          axis.ticks = element_blank())+
    scale_fill_manual(values=colors)+
    theme(axis.text.x = element_text(size=15,angle = 45, vjust = 1,hjust = 1),axis.text.y = element_text(size=15))+
    #geom_signif(comparisons = list(c("WT", "SCZ")),test.args = "less",map_signif_level=TRUE,size = 1,textsize = 8,color="black",vjust = 0.5 )+
    theme(legend.position="none")+
    annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 0.7)+
    annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
  return(plot)
}



