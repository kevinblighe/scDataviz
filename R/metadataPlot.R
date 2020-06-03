#' Colour shade a 2-dimensional reduction / embedding based on metadata, typically contained within a \code{SingleCellExperiment} object.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'   data-frame or matrix, columns named in \code{dimColnames} will be extracted
#'   from the data and used to generate the plot. If a SingleCellExperiment
#'   object, a reduction named by \code{reducedDim} will be taken from your object
#'   and used to generate the plot, again using columns whose names are
#'   specified in \code{dimColnames}.
#' @param meta If 'indata' is a non-SingleCellExperiment object, 'meta' must be
#'   activated and relate to a data-frame of metadata that aligns with the rows
#'   of \code{indata}, and that also contains a column name specified by \code{colby}.
#' @param reducedDim A reduced dimensional embedding stored within \code{indata},
#'   e.g., PCA or UMAP.
#' @param dimColnames The column names of the dimensions to use.
#' @param colby If NULL, all points will be coloured differently. If not NULL,
#'   the value is assumed to be a column name in \code{metadata(indata)} relating
#'   to some grouping / categorical variable.
#' @param colkey Vector of name-value pairs relating to value passed to 'col',
#'   e.g., \code{c(A='forestgreen', B='gold')}.
#' @param pointSize Size of plotted points.
#' @param legendPosition Position of legend \code{('top', 'bottom', 'left', 'right',
#'   'none')}.
#' @param legendLabSize Size of plot legend text.
#' @param legendIconSize Size of plot legend icons / symbols.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param celllab A vector containing any cells that the user wishes to label
#'   in the plot.
#' @param labSize Size of labels.
#' @param labhjust Horizontal adjustment of label.
#' @param labvjust Vertical adjustment of label.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding points by line connectors.
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
#'
#' @details
#' Colour shade a 2-dimensional reduction / embedding based on metadata, typically contained within a \code{SingleCellExperiment} object.
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
#' u <- umap::umap(mat)$layout
#' colnames(u) <- c('UMAP1','UMAP2')
#' rownames(u) <- rownames(mat)
#'
#' metadata <- data.frame(
#'   group = c(rep('PB1', 25), rep('PB2', 25)),
#'   row.names = rownames(u))
#'
#' metadataPlot(u, meta = metadata, colby = 'group')
#'
#' @import SingleCellExperiment ggplot2 ggrepel
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom methods is
#' 
#' @export
metadataPlot <- function(
  indata,
  meta = NULL,
  reducedDim = 'UMAP',
  dimColnames = c('UMAP1','UMAP2'),
  colby = NULL,
  colkey = NULL,
  pointSize = 0.5,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 5.0,
  xlim = NULL,
  ylim = NULL,
  celllab = NULL,
  labSize = 3.0,
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
  title = 'Metadata plot',
  subtitle = '',
  caption = ifelse(is(indata, 'SingleCellExperiment'),
    paste0('Total cells, ',
      nrow(as.data.frame(reducedDim(indata, reducedDim)))),
    paste0('Total cells, ', nrow(meta))),
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
  dim1 <- dim2 <- lab <- NULL

  # pull in the base theme, and add on parameters if necessary
  th <- basetheme(titleLabSize, subtitleLabSize, captionLabSize,
    axisLabSize, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, legendPosition, legendLabSize)

  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')
    plotobj <- as.data.frame(reducedDim(indata, reducedDim)[,dimColnames])
    plotobj <- data.frame(plotobj, metadata(indata))
    colnames(plotobj) <- c('dim1','dim2',colnames(metadata(indata)))

  } else {

    message('--input data class is ', class(indata))

    if (is.null(meta)) {
      stop('When the input data is a non-SingleCellExperiment object, ',
        '\'indata\' must relate to a 2-dimensional reduction / embedding that',
        'contains columns specified by \'dimColnames\', while \'meta\'',
        ' must be non-NULL and comprise a data-frame of metadata that ',
        'describes \'indata\', and which contains a column name',
        ' specified by \'colby\'.')
    } else if (!all(rownames(meta) == rownames(indata))) {
      stop('\'rownames(meta)\' must be equal to \'rownames(indata)\'')
    }

    plotobj <- as.data.frame(indata[,dimColnames])
    plotobj <- data.frame(plotobj, meta)
    colnames(plotobj) <- c('dim1','dim2', colnames(meta))
  }

  # set plot labels (e.g. cell names)
  if (!is.null(celllab)) {
    plotobj$lab <- rownames(plotobj)
    plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)

    names.new <- rep(NA, length(plotobj$lab))
    indices <- which(plotobj$lab %in% celllab)
    names.new[indices] <- plotobj$lab[indices]
    plotobj$lab <- names.new
  }

  plotobj$col <- plotobj[,colby]

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

  # create the plot object
  plot <- ggplot(plotobj, aes(x = dim1, y = dim2)) + th +

    guides(fill = guide_legend(),
      colour = guide_legend(override.aes = list(size = legendIconSize)))

  plot <- plot + geom_point(aes(color = col),
        size = pointSize)

  # sort out custom colour pairing,
  if (!is.null(colkey)) {
    plot <- plot + scale_colour_discrete('') +
      scale_color_manual(values = colkey)
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

  if (!is.null(celllab)) {
    if (drawConnectors && is.null(celllab)) {
      plot <- plot + geom_text_repel(
        data = plotobj,
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (drawConnectors && !is.null(celllab)) {
      plot <- plot + geom_text_repel(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (!drawConnectors && !is.null(celllab)) {
      plot <- plot + geom_text(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
    } else if (!drawConnectors && is.null(celllab)) {
      plot <- plot + geom_text(
        data = plotobj,
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
    }
  }

  return(plot)
}
