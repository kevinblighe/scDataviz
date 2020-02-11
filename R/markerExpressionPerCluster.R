markerExpressionPerCluster <- function(
  indata,
  assay = 'scaled',
  clusters = sample(unique(metadata(indata)[['Cluster']]), 9),
  clusterAssign = metadata(indata)[['Cluster']],
  markers = sample(rownames(indata), 10),
  ncol = 5,
  nrow = 2,
  legendPosition = 'none',
  legendLabSize = 12,
  legendIconSize = 5.0,
  legendKeyHeight = 2.5,
  xlim = NULL,
  ylim = NULL,
  yfixed = FALSE,
  xlab = 'Marker',
  xlabAngle = 90,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = 'Expression',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  stripLabSize = 16,
  title = 'Marker expression per cluster',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  borderWidth = 0.8,
  borderColour = 'black')
{
  Marker <- Expression <- lab <- NULL

  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background = element_rect(),

      title = element_text(size = legendLabSize),

      plot.title=element_text(angle=0, size=titleLabSize,
        face='bold', vjust=1),
      plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
        face = 'plain', vjust = 1),
      plot.caption=element_text(angle = 0, size = captionLabSize,
        face = 'plain', vjust = 1),

      axis.text.x=element_text(angle = xlabAngle, size = axisLabSize,
        hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = axisLabSize,
        hjust = ylabhjust, vjust = ylabvjust),
      axis.title = element_text(size = axisLabSize),

      legend.title = element_blank(),
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = legendLabSize),
      legend.key.height = unit(legendKeyHeight, 'cm'),

      strip.text.x = element_text(size = stripLabSize, face = 'bold', margin = margin(b = 5, t = 5)))

  if (class(indata) == 'SingleCellExperiment') {

    message('--input data class is SingleCellExperiment')
    plotobj <- data.frame(Cluster = clusterAssign,
      as.data.frame(t(as.matrix(assay(indata, assay)[which(rownames(indata) %in% markers),]))))
    plotobj <- plotobj[which(plotobj$Cluster %in% clusters),]
    plotobj <- melt(plotobj, id.vars = 'Cluster')
    colnames(plotobj) <- c('Cluster','Marker','Expression')
    plotobj$Cluster <- paste0('Cluster ', plotobj$Cluster)
    plotobj$Cluster <- factor(plotobj$Cluster, levels = unique(paste0('Cluster ', sort(clusters))))

  } else {

    message('--input data class is ', class(indata))
    plotobj <- data.frame(Cluster = clusterAssign,
      as.data.frame(t(as.matrix(indata[which(rownames(indata) %in% markers),]))))
    plotobj <- plotobj[which(plotobj$Cluster %in% clusters),]
    plotobj <- melt(plotobj, id.vars = 'Cluster')
    colnames(plotobj) <- c('Cluster','Marker','Expression')
    plotobj$Cluster <- factor(plotobj$Cluster, levels = unique(sort(clusters)))

  }

  # initialise the plot object
  plot <- ggplot(plotobj, aes(x = Marker, y = Expression)) + th +

    guides(
      fill = guide_legend(),
      shape = guide_legend(),
      alpha = FALSE)

  plot <- plot + geom_boxplot(
    position = position_dodge(width = 0.1),
    outlier.shape = '.',
    outlier.colour = 'black',
    outlier.size = 0.1,
    aes(fill = Marker))

  if (yfixed == TRUE) {
    plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol)
  } else {
    plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol, scales = 'free_y')
  }

  # add elements to the plot for xy labeling and axis limits
  plot <- plot + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # border around plot
  plot <- plot +
    theme(panel.border = element_rect(
      colour = borderColour,
      fill = NA,
      size = borderWidth))

  return(plot)
}
