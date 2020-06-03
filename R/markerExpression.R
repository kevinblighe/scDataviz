#' Highlight the individual marker expression profile across a 2-dimensional reduction / embedding, typically contained within a \code{SingleCellExperiment} object. By default, this function plots the expression profile of 6 randomly-selected markers from your data.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'   data-frame or matrix, this should relate to expression data (cells as
#'   columns; genes as rows). If a \code{SingleCellExperiment} object, data will be
#'   extracted from an assay component named by \code{assay}.
#' @param layout If 'indata' is a non-SingleCellExperiment object, \code{layout} must
#'   be activated and relate to a 2-dimensional reduction / embedding, although,
#'   technically, any data-frame or matrix of numbers will be accepted, provided
#'   that it aligns with the dimensions of \code{indata}, and provided that it
#'   contains columns as specified in \code{dimColnames}.
#' @param assay Name of the assay slot in 'indata' from which data will be
#'   taken, assuming \code{indata} is a \code{SingleCellExperiment} object.
#' @param reducedDim A reduced dimensional component stored within \code{indata},
#'   e.g., PCA or UMAP.
#' @param dimColnames The column names of the dimensions to use.
#' @param markers Vector containing marker names to plot.
#' @param ncol Number of columns for faceting.
#' @param nrow Number of rows for faceting.
#' @param col Colours used for generation of fill gradient according to
#'   expression values. Can be 2 or 3 colours.
#' @param colMidpoint Mid-point (expression value) for the colour range. Only
#'   used when 3 colours are specified by \code{col}.
#' @param alpha Control the gradient of colour transparency, with 1 being opaque.
#' @param pointSize Size of plotted points.
#' @param legendPosition Position of legend \code{('top', 'bottom', 'left', 'right',
#'   'none')}.
#' @param legendLabSize Size of plot legend text.
#' @param legendKeyHeight Height of the legend key.
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
#' @param stripLabSize Size of the strip (marker) labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param hline Draw one or more horizontal lines passing through this/these
#'    values on y-axis. For single values, only a single numerical value is
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
#' Highlight the individual marker expression profile across a 2-dimensional reduction / embedding, typically contained within a \code{SingleCellExperiment} object. By default, this function plots the expression profile of 6 randomly-selected markers from your data.
#'
#' @return A \code{ggplot2} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#' MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#' rownames(mat) <- paste0('cell', 1:nrow(mat))
#'
#' u <- umap::umap(mat)$layout
#' colnames(u) <- c('UMAP1','UMAP2')
#' rownames(u) <- rownames(mat)
#'
#' markerExpression(t(mat), layout = u)
#'
#' @import SingleCellExperiment ggplot2
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom reshape2 melt
#' @importFrom methods is
#' 
#' @export
markerExpression <- function(
  indata,
  layout = NULL,
  assay = 'scaled',
  reducedDim = 'UMAP',
  dimColnames = c('UMAP1','UMAP2'),
  markers = sample(rownames(indata), 6),
  ncol = 3,
  nrow = 2,
  col = c('darkblue', 'yellow'),
  colMidpoint = 0,
  alpha = c(0.0, 1),
  pointSize = 0.5,
  legendPosition = 'right',
  legendLabSize = 12,
  legendKeyHeight = 2.5,
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
  stripLabSize = 16,
  title = 'Individual marker expression',
  subtitle = '',
  caption = ifelse(is(indata, 'SingleCellExperiment'),
    paste0('Total cells, ',
      nrow(as.data.frame(reducedDim(indata, reducedDim)))),
    paste0('Total cells, ', nrow(layout))),
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
  dim1 <- dim2 <- Expression <- lab <- NULL

  # pull in the base theme, and add on parameters if necessary
  th <- basetheme(titleLabSize, subtitleLabSize, captionLabSize,
    axisLabSize, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, legendPosition, legendLabSize) +

    theme(legend.key.height = unit(legendKeyHeight, 'cm'),
      strip.text.x = element_text(size = stripLabSize,
        face = 'bold', margin = margin(b = 5, t = 5)))

  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')
    plotobj <- as.data.frame(reducedDim(indata, reducedDim)[,dimColnames])
    plotobj <- data.frame(plotobj,
      as.data.frame(t(as.matrix(assay(indata, assay)))))
    plotobj <- melt(plotobj, id.vars = dimColnames)

  } else {

    message('--input data class is ', class(indata))

    if (is.null(layout)) {
      stop('When the input data is a non-SingleCellExperiment object, ',
        '\'indata\' must relate to an expression matrix (cells as columns; ',
        'genes as rows), while \'layout\' must be non-NULL and relate to ',
        'a 2-dimensional embedding containing columns specified by ',
        '\'dimColnames\'')
    } else if (!all(rownames(layout) == colnames(indata))) {
      stop('\'rownames(layout)\' must be equal to \'colnames(indata)\'')
    }

    plotobj <- as.data.frame(layout[,dimColnames])
    plotobj <- data.frame(plotobj, as.data.frame(t(as.matrix(indata))))
    plotobj <- melt(plotobj, id.vars = dimColnames)
  }
  colnames(plotobj) <- c('dim1','dim2','Marker','Expression')

  plotobj <- plotobj[which(plotobj$Marker %in% markers),]

  # set plot labels (e.g. cell names)
  if (!is.null(celllab)) {
    plotobj$lab <- rep(colnames(indata), length(markers))
    plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)

    names.new <- rep(NA, length(plotobj$lab))
    indices <- which(plotobj$lab %in% celllab)
    names.new[indices] <- plotobj$lab[indices]
    plotobj$lab <- names.new
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

  # order by expression level to ensure that highly expressed are coloured last
  plotobj <- plotobj[order(plotobj$Expression, decreasing = FALSE),]

  # initialise the plot object
  plot <- ggplot(plotobj, aes(x = dim1, y = dim2, alpha = Expression)) + th +

    guides(
      fill = guide_legend(),
      shape = guide_legend(),
      alpha = FALSE)

  plot <- plot + geom_point(aes(colour = Expression), size = pointSize)

  if (length(col) == 2) {
    plot <- plot +
       scale_colour_gradient(
         low = col[1],
         high = col[2],
         name = 'Expression')
  } else if (length(col) == 3) {
    plot <- plot +
      scale_colour_gradient2(
        low = col[1],
        mid = col[2],
        high = col[3],
        midpoint = colMidpoint,
        limits = c(min(plotobj$Expression), max(plotobj$Expression)),
        space = 'Lab',
        name = 'Expression')
  }

  plot <- plot + #scale_alpha(range = c(alpha[1], alpha[2]), guide = FALSE) +
    facet_wrap( ~ Marker, nrow = nrow, ncol = ncol)

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
