\name{markerExpressionPerCluster}

\alias{markerExpressionPerCluster}

\title{markerExpressionPerCluster}

\description{Generate box-and-whisker plots illustrating marker expression per k-NN identified cluster. By default, 9 randomly-selected clusters are selected, and the expression profiles of 10 randomly-selected markers are plot across these.}

\usage{
  markerExpressionPerCluster(
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
    borderColour = 'black')
}

\arguments{
  \item{indata}{A data-frame or matrix, or SingleCellExperiment object. If a
    data-frame or matrix, this should relate to expression data (cells as
    columns; genes as rows). If a SingleCellExperiment object, data will be
    extracted from an assay component named by 'assay'. REQUIRED.}
  \item{assay}{Name of the assay slot in 'indata' from which data will be
    taken, assuming 'indata' is a SingleCellExperiment object.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{clusters}{Vector containing clusters to plot. DEFAULT =
    sample(unique(metadata(indata)[['Cluster']]), 5). OPTIONAL.}
  \item{clusterAssign}{A vector of cell-to-cluster assignments. This can be
    from any source but must align with your cells / variables. There is no
    check to ensure this when 'indata' is not a SingleCellExperiment object.
    DEFAULT = metadata(indata)[['Cluster']]. OPTIONAL.}
  \item{markers}{Vector containing marker names to plot.
    Default = sample(rownames(indata), 10). OPTIONAL.}
  \item{ncol}{Number of columns for faceting. DEFAULT = 5. OPTIONAL.}
  \item{nrow}{Number of rows for faceting. DEFAULT = 2. OPTIONAL.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
  'none'). DEFAULT = 'none'. OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 12. OPTIONAL.}
  \item{legendKeyHeight}{Height of the legend key. DEFAULT = 2.5. OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = NULL. OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = NULL. OPTIONAL.}
  \item{yfixed}{Logical, specifying whether or not to fix the y-axis
    scales across all clusters when faceting. DEFAULT = FALSE. OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = 'Marker'. OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 90. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5. OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = 'Expression'. OPTIONAL.}
  \item{ylabAngle}{Rotation angle of y-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{ylabhjust}{Horizontal adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylabvjust}{Vertical adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{axisLabSize}{Size of x- and y-axis labels. DEFAULT = 16. OPTIONAL.}
  \item{stripLabSize}{Size of the strip labels. DEFAULT = 16. OPTIONAL.}
  \item{title}{Plot title. DEFAULT = 'Marker expression per cluster'. OPTIONAL.}
  \item{subtitle}{Plot subtitle. DEFAULT = ''. OPTIONAL.}
  \item{caption}{Plot caption. DEFAULT = ''. OPTIONAL.}
  \item{titleLabSize}{Size of plot title. DEFAULT = 16. OPTIONAL.}
  \item{subtitleLabSize}{Size of plot subtitle. DEFAULT = 12. OPTIONAL.}
  \item{captionLabSize}{Size of plot caption. DEFAULT = 12. OPTIONAL.}
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
  # create random data that follows a negative binomial
  mat <- jitter(matrix(
    MASS::rnegbin(rexp(5000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat) <- paste0('CD', 1:ncol(mat))
  rownames(mat) <- paste0('cell', 1:nrow(mat))

  clus <- clusKNN(mat)
  markerExpressionPerCluster(t(mat), clusters = c(0, 1), clusterAssign = clus)
}
