#' Highlight cell-to-cluster assignments across a 2-dimensional reduction / embedding.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'   data-frame or matrix, columns named in \code{dimColnames} will be extracted
#'   from the data and used to generate the plot. If a SingleCellExperiment
#'   object, a reduction named by \code{reducedDim} will be taken from your object
#'   and used to generate the plot, again using columns whose names are
#'   specified in \code{dimColnames}.
#' @param clusterVector If \code{indata} is a non-\code{SingleCellExperiment} object,
#'   \code{clusterVector} must be non-NULL and relate to a cell-to-cluster
#'   assignment whose length matches \code{nrow(indata)}.
#' @param reducedDim A reduced dimensional embedding stored within 'indata',
#'   e.g., PCA or UMAP.
#' @param dimColnames The column names of the dimensions to use.
#' @param clusterColname The column name in the metadata of \code{indata} that
#'   contains the cell-to-cluster assignment, assuming \code{indata} is a
#'   \code{SingleCellExperiment} object.
#' @param pointSize Size of plotted points.
#' @param legendPosition Position of legend \code{('top', 'bottom', 'left', 'right',
#'   'none')}.
#' @param legendLabSize Size of plot legend text.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param label Logical, indicating whether or not to label the clusters.
#' @param labSize Size of labels.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding cluster islands by line connectors.
#' @param widthConnectors Line width of connectors.
#' @param colConnectors Line colour of connectors.
#' @param xlab Label for x-axis.
#' @param xlabAngle Rotation angle of x-axis labels.
#' @param xlabhjust Horizontal adjustment of x-axis labels.
#' @param xlabvjust Vertical adjustment of x-axis labels.
#' @param ylab Label for y-axis.
#' @param ylabAngle Rotation angle of y-axis labels.
#' @param ylabhjust Horizontal adjustment of y-axis labels.
#' @param ylabvjust Vertical adjustment of y-axis labels.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param hline Draw one or more horizontal lines passing through this/these
#'   values on y-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param hlineType Line type for hline \code{('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash')}.
#' @param hlineCol Colour of hline.
#' @param hlineWidth Width of hline.
#' @param vline Draw one or more vertical lines passing through this/these
#'   values on x-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param vlineType Line type for vline \code{('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash')}.
#' @param vlineCol Colour of vline.
#' @param vlineWidth Width of vline.
#' @param gridlines.major Logical, indicating whether or not to draw major
#'   gridlines.
#' @param gridlines.minor Logical, indicating whether or not to draw minor
#'   gridlines.
#' @param borderWidth Width of the border on the x and y axes.
#' @param borderColour Colour of the border on the x and y axes.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Highlight cell-to-cluster assignments across a 2-dimensional reduction / embedding.
#'
#' @return A \code{ggplot2} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#' rownames(mat) <- paste0('cell', 1:nrow(mat))
#'
#' u <- umap::umap(mat)
#' clusvec <- clusKNN(u$layout)
#' plotClusters(u$layout, clusvec)
#'
#' @import SingleCellExperiment ggplot2 ggrepel
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom methods is
#'
#' @export
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
  borderColour = 'black',
  verbose = TRUE)
{
  dim1 <- dim2 <- Cluster <- lab <- NULL

  # pull in the base theme, and add on parameters if necessary
  th <- basetheme(titleLabSize, subtitleLabSize, captionLabSize,
    axisLabSize, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, legendPosition, legendLabSize)

  if (is(indata, 'SingleCellExperiment')) {

    if (verbose) message('--input data class is SingleCellExperiment')
    plotobj <- as.data.frame(reducedDim(indata, reducedDim))
    plotobj <- data.frame(plotobj,
      Cluster = metadata(indata)[[clusterColname]])
    colnames(plotobj) <- c('dim1','dim2','Cluster')
    plotobj$Cluster <- factor(plotobj$Cluster)

  } else {

    if (verbose) message('--input data class is ', class(indata))

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
          xlim = c(NA, NA),
          ylim = c(NA, NA),
          direction = 'both',
          max.overlaps = Inf,
          min.segment.length = 0,
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          na.rm = TRUE)
    } else if (!drawConnectors) {
      plot <- plot + geom_text(
        data = subset(plotobj, !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE)
    }
  }

  return(plot)
}
