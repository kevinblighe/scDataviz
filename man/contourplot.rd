\name{contourplot}

\alias{contourplot}

\title{contourplot}

\description{Draw a contour plot representation of a 2-dimensional embedding.}

\usage{
  contourplot(sce,
  reducedDim = 'UMAP',
  lowcol = 'darkblue',
  highcol = 'darkred',
  alpha = c(0.0, 0.5),
  contour = 'black',
  bins = 300,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 5.0,
  legendKeyHeight = 2.5,
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
  title = 'Cellular density and contours',
  subtitle = '',
  caption = paste0('Total cells, ', nrow(as.data.frame(reducedDim(sce, "UMAP"))), '; Bins, ', bins),
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
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.},
  \item{reducedDim}{A reduced dimensional component stored within 'sce',
    e.g., PCA or UMAP. DEFAULT = 'UMAP'. OPTIONAL.}
  \item{lowcol}{Shade for low-density contours. DEFAULT = 'darkblue'. OPTIONAL.}
  \item{highcol}{Shade for high-density contours. DEFAULT = 'darkred'. OPTIONAL.}
  \item{alpha}{Control the gradient of colour transparency, with 1 being opaque.
    DEFAULT = c(0.0, 0.5). OPTIONAL.}
  \item{contour}{The colour of the contour lines. DEFAULT = 'black'. OPTIONAL.}
  \item{bins}{The number of bins that determine the overall density values.
    DEFAULT = 300. OPTIONAL.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
  'none'). DEFAULT = 'right'. OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 12. OPTIONAL.}
  \item{legendIconSize}{Size of plot legend icons / symbols. DEFAULT = 5.0.
  OPTIONAL.}
  \item{legendKeyHeight}{Height of the legend key. DEFAULT = 2.5. OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = NULL. OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = NULL. OPTIONAL.}
  \item{celllab}{A vector containing any cells that the user wishes to label
    in the plot. DEFAULT = NULL. OPTIONAL.}
  \item{labSize}{Size of labels. DEFAULT = 3.0. OPTIONAL.}
  \item{labhjust}{Horizontal adjustment of label. DEFAULT = 1.5. OPTIONAL.}
  \item{labvjust}{Vertical adjustment of label. DEFAULT = 0. OPTIONAL.}

  \item{drawConnectors}{Logical, indicating whether or not to connect plot
  labels to their corresponding points by line connectors. DEFAULT = TRUE.
  OPTIONAL.}
  \item{widthConnectors}{Line width of connectors. DEFAULT = 0.5. OPTIONAL.}
  \item{colConnectors}{Line colour of connectors. DEFAULT = 'grey50'. OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = 'UMAP1'. OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5. OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = 'UMAP2'. OPTIONAL.}
  \item{ylabAngle}{Rotation angle of y-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{ylabhjust}{Horizontal adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylabvjust}{Vertical adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{axisLabSize}{Size of x- and y-axis labels. DEFAULT = 16. OPTIONAL.}
  \item{title}{Plot title. DEFAULT = 'Cellular density and contours'. OPTIONAL.}
  \item{subtitle}{Plot subtitle. DEFAULT = ''. OPTIONAL.}
  \item{caption}{Plot caption. DEFAULT =
    paste0('Total cells, ', nrow(as.data.frame(reducedDim(sce, "UMAP"))),
    '; Bins, ', bins). OPTIONAL.}
  \item{titleLabSize}{Size of plot title. DEFAULT = 16. OPTIONAL.}
  \item{subtitleLabSize}{Size of plot subtitle. DEFAULT = 12. OPTIONAL.}
  \item{captionLabSize}{Size of plot caption. DEFAULT = 12. OPTIONAL.}
  \item{hline}{Draw one or more horizontal lines passing through this/these
  values on y-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = NULL. OPTIONAL.}
  \item{hlineType}{Line type for hline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{hlineCol}{Colour of hline. DEFAULT = 'black'. OPTIONAL.}
  \item{hlineWidth}{Width of hline. DEFAULT = 0.4. OPTIONAL.}
  \item{vline}{Draw one or more vertical lines passing through this/these
  values on x-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = NULL. OPTIONAL.}
  \item{vlineType}{Line type for vline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{vlineCol}{Colour of vline. DEFAULT = 'black'. OPTIONAL.}
  \item{vlineWidth}{Width of vline. DEFAULT = 0.4. OPTIONAL.}
  \item{gridlines.major}{Logical, indicating whether or not to draw major
  gridlines. DEFAULT = TRUE. OPTIONAL.}
  \item{gridlines.minor}{Logical, indicating whether or not to draw minor
  gridlines. DEFAULT = TRUE. OPTIONAL.}
  \item{borderWidth}{Width of the border on the x and y axes. DEFAULT = 0.8.
  OPTIONAL.}
  \item{borderColour}{Colour of the border on the x and y axes. DEFAULT =
  'black'. OPTIONAL.}
}

\value{
A \code{\link{ggplot2}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  # not run
}
