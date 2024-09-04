plotCorrelation <-function(x, y, groups=NULL,legend.pos = "none",maintitle=NULL) {
  
  pval <- paste(sprintf("Correlation = %.3f\nP", Hmisc::rcorr(x, y)$r[2,1]),
                paste0("= ", signif(Hmisc::rcorr(x, y)$P[2,1], 3)) )
  cor_plot <- data.frame(a=x, b=y)
  plot <- ggplot(cor_plot, aes(x=a, y=b)) +
    geom_smooth(method='lm',formula=y~x, se=F, linetype = "dashed", colour="gray") +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = pval)+
    geom_point(aes(color=NULL)) +ggtitle(maintitle)
  plot
}