#' Generate box-and-whisker plots illustrating marker expression per k-NN identified cluster. By default, 5 randomly-selected clusters are selected, and the expression profiles of 10 randomly-selected markers are plot across these.
#'
#' @param indata A data-frame or matrix, or SingleCellExperiment object. If a
#'    data-frame or matrix, this should relate to expression data (cells as
#'    columns; genes as rows). If a \code{SingleCellExperiment} object, data will be
#'    extracted from an assay component named by \code{assay}.
#' @param assay Name of the assay slot in \code{indata} from which data will be
#'    taken, assuming \code{indata} is a \code{SingleCellExperiment} object.
#' @param clusters Vector containing clusters to plot.
#' @param clusterAssign A vector of cell-to-cluster assignments. This can be
#'    from any source but must align with your cells / variables. There is no
#'    check to ensure this when \code{indata} is not a \code{SingleCellExperiment} object.
#' @param markers Vector containing marker names to plot.
#' @param ncol Number of columns for faceting.
#' @param nrow Number of rows for faceting.
#' @param legendPosition Position of legend \code{('top', 'bottom', 'left', 'right',
#'  'none')}.
#' @param legendLabSize Size of plot legend text.
#' @param legendKeyHeight Height of the legend key.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param yfixed Logical, specifying whether or not to fix the y-axis
#'    scales across all clusters when faceting.
#' @param xlab Label for x-axis.
#' @param xlabAngle Rotation angle of x-axis labels.
#' @param xlabhjust Horizontal adjustment of x-axis labels.
#' @param xlabvjust Vertical adjustment of x-axis labels.
#' @param ylab Label for y-axis.
#' @param ylabAngle Rotation angle of y-axis labels.
#' @param ylabhjust Horizontal adjustment of y-axis labels.
#' @param ylabvjust Vertical adjustment of y-axis labels.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param stripLabSize Size of the strip labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param borderWidth Width of the border on the x and y axes.
#' @param borderColour Colour of the border on the x and y axes.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Generate box-and-whisker plots illustrating marker expression per k-NN identified cluster. By default, 5 randomly-selected clusters are selected, and the expression profiles of 10 randomly-selected markers are plot across these.
#'
#' @return A \code{ggplot2} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(5000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#' rownames(mat) <- paste0('cell', 1:nrow(mat))
#'
#' clus <- clusKNN(mat)
#' markerExpressionPerCluster(t(mat), clusters = c(0, 1),
#'   clusterAssign = clus)
#'
#' @import SingleCellExperiment ggplot2
#'
#' @importFrom MASS rnegbin
#' 
#' @export
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
