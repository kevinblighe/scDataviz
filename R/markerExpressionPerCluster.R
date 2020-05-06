markerExpressionPerCluster <- function(
  indata,
  assay = 'scaled',
  clusters = sample(unique(metadata(indata)[['Cluster']]), 5),
  clusterAssign = metadata(indata)[['Cluster']],
  markers = sample(rownames(indata), 10),
  ncol = 5,
  nrow = 2,
  legendPosition = 'none',
  legendLabSize = 12,
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
  borderColour = 'black',
  verbose = TRUE)
{
  Marker <- Expression <- lab <- NULL

  # pull in the base theme, and add on parameters if necessary
  th <- basetheme(titleLabSize, subtitleLabSize, captionLabSize,
    axisLabSize, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, legendPosition, legendLabSize) +

    theme(legend.key.height = unit(legendKeyHeight, 'cm'),
      strip.text.x = element_text(size = stripLabSize,
        face = 'bold', margin = margin(b = 5, t = 5)))

  if (is(indata, 'SingleCellExperiment')) {

    if (verbose) message('--input data class is SingleCellExperiment')
    idx <- which(rownames(indata) %in% markers)
    plotobj <- data.frame(Cluster = clusterAssign,
      as.data.frame(t(as.matrix(assay(indata, assay)[idx,]))))
    plotobj <- plotobj[which(plotobj$Cluster %in% clusters),]
    plotobj <- melt(plotobj, id.vars = 'Cluster')
    colnames(plotobj) <- c('Cluster','Marker','Expression')
    plotobj$Cluster <- paste0('Cluster ', plotobj$Cluster)
    plotobj$Cluster <- factor(plotobj$Cluster,
      levels = unique(paste0('Cluster ', sort(clusters))))

  } else {

    if (verbose) message('--input data class is ', class(indata))
    idx <- which(rownames(indata) %in% markers)
    plotobj <- data.frame(Cluster = clusterAssign,
      as.data.frame(t(as.matrix(indata[idx,]))))
    plotobj <- plotobj[which(plotobj$Cluster %in% clusters),]
    plotobj <- melt(plotobj, id.vars = 'Cluster')
    colnames(plotobj) <- c('Cluster','Marker','Expression')
    plotobj$Cluster <- factor(plotobj$Cluster,
      levels = unique(sort(clusters)))

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

  if (yfixed) {
    plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol)
  } else {
    plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol,
      scales = 'free_y')
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
