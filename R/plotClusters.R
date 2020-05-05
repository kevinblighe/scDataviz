plotClusters <- function(
  indata,
  clusterVector = NULL,
  reducedDim = 'UMAP',
  dimColnames = c('UMAP1','UMAP2'),
  clusterColname = 'Cluster',
  pointSize = 0.5,
  legendPosition = 'none',
  legendLabSize = 12,
  xlim = NULL,
  ylim = NULL,
  label = TRUE,
  labSize = 5.0,
  labhjust = 1.5,
  labvjust = 0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  xlab = dimColnames[1],
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = dimColnames[2],
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  title = 'k-nearest neighbour (k-NN) clusters',
  subtitle = '',
  caption = ifelse(is(indata, 'SingleCellExperiment'),
    paste0('Total cells, ',
      nrow(as.data.frame(reducedDim(indata, reducedDim)))),
    paste0('Total cells, ', length(clusterVector))),
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  borderWidth = 0.8,
  borderColour = 'black')
{
  metadata <- dim1 <- dim2 <- Cluster <- lab <- NULL

  # pull in the base theme, and add on parameters if necessary
  th <- basetheme(titleLabSize, subtitleLabSize, captionLabSize,
    axisLabSize, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, legendPosition, legendLabSize)

  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')
    plotobj <- as.data.frame(reducedDim(indata, reducedDim))
    plotobj <- data.frame(plotobj,
      Cluster = metadata(indata)[[clusterColname]])
    colnames(plotobj) <- c('dim1','dim2','Cluster')
    plotobj$Cluster <- factor(plotobj$Cluster)

  } else {

    message('--input data class is ', class(indata))

    if (is.null(clusterVector)) {
      stop('When the input data is a non-SingleCellExperiment object, ',
        '\'indata\' must relate to 2-dimensional embedding, while ',
        '\'clusterVector\' must be non-NULL and relate to a vector ',
        'of cell-to-cluster assignments whose length matches ',
        '\'nrow(indata)\'')
    }

    plotobj <- as.data.frame(indata)
    plotobj <- data.frame(plotobj, Cluster = clusterVector)
    colnames(plotobj) <- c('dim1','dim2','Cluster')
    plotobj$Cluster <- factor(plotobj$Cluster)

  }

  # set labels
  if (label == TRUE) {
    plotobj$lab <- paste0('Cluster', plotobj$Cluster)
    plotobj$lab[duplicated(plotobj$lab)] <- NA
  }

  if (is.null(xlim)) {
    xlim <- c(
      min(plotobj[,'dim1'], na.rm = TRUE) - 1,
      max(plotobj[,'dim1'], na.rm = TRUE) + 1)
  }

  if (is.null(ylim)) {
    ylim <- c(
      min(plotobj[,'dim2'], na.rm = TRUE) - 1,
      max(plotobj[,'dim2'], na.rm = TRUE) + 1)
  }

  # initialise the plot object
  plot <- ggplot(plotobj, aes(x = dim1, y = dim2)) + th +

    guides(
      fill = guide_legend(),
      shape = guide_legend(),
      alpha = FALSE)

  plot <- plot + geom_point(aes(colour = Cluster), size = pointSize)

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

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = hline,
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # border around plot
  plot <- plot +
    theme(panel.border = element_rect(
      colour = borderColour,
      fill = NA,
      size = borderWidth))

  # gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  if (label) {
    if (drawConnectors) {
      plot <- plot + geom_text_repel(
        data = subset(plotobj, !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (!drawConnectors) {
      plot <- plot + geom_text(
        data = subset(plotobj, !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
    }
  }

  return(plot)
}
