markerExpressionPerCluster <- function(
  data,
  cluster = c(0:length(unique(data$nnc[,1]))),
  ncol = 2,
  markernames = NULL)
{
  mat <- data$expression
  ggdata <- melt(data.frame(data$nnc, mat))
  colnames(ggdata) <- c('Cluster', 'Marker', 'Expression')
  ggdata$Marker <- factor(ggdata$Marker, levels = make.names(colnames(mat)))

  if (!is.null(markernames)) {
    levels(ggdata$Marker) <- markernames
  }

  ggdata$Cluster <- factor(ggdata$Cluster, levels = cluster)
  ggdata <- ggdata[!is.na(ggdata$Cluster),]

  gg <- ggplot(data = ggdata, aes(x = Marker, y = Expression)) +

    geom_boxplot(position=position_dodge(width=0.1), outlier.shape = ".", outlier.colour="black", outlier.size = 0.1, aes(fill = Marker)) +

    facet_wrap(~ Cluster, ncol = ncol, scales = 'free_y') +

    #Set the size of the plotting window
    theme_bw(base_size=24) +

    #Modify various aspects of the plot text and legend
    theme(
      legend.position="none",
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=18, face="bold", vjust=1),

      axis.text.x=element_text(angle=45, size=12, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
      axis.title=element_text(size=18, face="bold"),

      #Legend
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=12),  #Text size
      title=element_text(size=12),      #Title text size

      strip.text.x = element_text(size=18, face="bold", margin = margin(b = 0, t = 0))) +

    #Change the size of the icons/symbols in the legend
    guides(colour=guide_legend(override.aes=list(size=2.5))) +

    #Set x- and y-axes labels
    xlab("") +
    ylab("Expression") +

    ggtitle('Cluster marker expression')

  return(gg)
}