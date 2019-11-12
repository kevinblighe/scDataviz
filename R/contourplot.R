contourplot <- function(
  sce,
  lowcol = 'darkblue',
  highcol = 'darkred',
  alpha = c(0.0, 0.5),
  contour = 'black',
  bins = 300,

  legendPosition = 'none',
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
  colConnectors = 'grey50',

  xlab = 'UMAP1',
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,

  ylab = 'UMAP2',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,

  title = NULL,
  subtitle = NULL,
  caption = NULL,
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
  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background=element_rect(),

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
      axis.title=element_text(size=axisLabSize),

      legend.position=legendPosition,
      legend.key=element_blank(),
      legend.key.size=unit(0.5, 'cm'),
      legend.text=element_text(size=legendLabSize),

      title=element_text(size=legendLabSize),
      legend.title=element_blank())

  plotobj <- as.data.frame(reducedDim(sce, "UMAP"))

  # set plot labels (e.g. cell names)
  if (!is.null(celllab)) {
    plotobj$lab <- rownames(plotobj)
    plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)

    names.new <- rep(NA, length(plotobj$lab))
    indices <- which(plotobj$lab %in% celllab)
    names.new[indices] <- plotobj$lab[indices]
    plotobj$lab <- names.new
  }

  if (is.null(title)) {
    title <- 'Cellular density and contours'
  }

  if (is.null(caption)) {
    caption <- paste0('Total cells, ', nrow(data), '; Bins, ', bins)
  }

  if (is.null(xlim)) {
    xlim <- c(
      min(plotobj[,'UMAP1'], na.rm = TRUE) - 1,
      max(plotobj[,'UMAP1'], na.rm = TRUE) + 1)
  }

  if (is.null(ylim)) {
    ylim <- c(
      min(plotobj[,'UMAP2'], na.rm = TRUE) - 1,
      max(plotobj[,'UMAP2'], na.rm = TRUE) + 1)
  }

  # initialise the plot object
  plot <- ggplot(plotobj, aes(UMAP1, UMAP2)) + th +

    guides(fill = guide_legend(),
      shape = guide_legend(),
      colour = guide_legend(override.aes = list(size = legendIconSize))) +

    stat_density2d(aes(alpha = ..level.., fill = ..level..), size = 1, bins = bins, geom = 'polygon') +

    scale_fill_gradient(low = lowcol, high = highcol) +

    scale_alpha(range = c(alpha[1], alpha[2]), guide = FALSE) +

    geom_density2d(colour = contour)

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
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }


  if (!is.null(celllab)) {
    if (drawConnectors == TRUE && is.null(celllab)) {
      plot <- plot + geom_text_repel(
        data = plotobj,
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (drawConnectors == TRUE && !is.null(celllab)) {
      plot <- plot + geom_text_repel(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (drawConnectors == FALSE && !is.null(celllab)) {
      plot <- plot + geom_text(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          check_overlap = TRUE,
          hjust = labhjust,
          vjust = labvjust)
    } else if (drawConnectors == FALSE && is.null(celllab)) {
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
