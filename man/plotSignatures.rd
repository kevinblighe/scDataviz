\name{plotSignatures}

\alias{plotSignatures}

\title{plotSignatures}

\description{Find enriched markers per identified cluster and visualise these as a custom corrplot.}

\usage{
  plotSignatures(sce,
  assay = 'scaled',
  clusterVector = metadata(sce)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  col = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100),
  cexlab = 1.0,
  cexlegend = 1.0,
  labDegree = 90)
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.}
  \item{assay}{Name of the assay slot in sce from which data will be taken.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{clusterVector}{A vector of cell-to-cluster assignments. This can be
    from any source but ought to be taken from the metadata. DEFAULT =
    metadata(sce)[['Cluster']]. OPTIONAL.}
  \item{funcSummarise}{A mathematical function used to summarise expression
    per marker per cluster. DEFAULT = function(x) median(x, na.rm = TRUE).
    OPTIONAL.}
  \item{col}{colorRampPalette to be used for shading low-to-high expression.
    DEFAULT = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100). OPTIONAL.}
  \item{cexlab}{cex of the main plot labels. DEFAULT = 1.0. OPTIONAL.}
  \item{cexlegend}{cex of the legend labels. DEFAULT = 1.0. OPTIONAL.}
  \item{labDegree}{Rotation angle of the main plot labels. DEFAULT = 90.
    OPTIONAL.}
}

\value{
A \code{\link{corrplot}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  # create random data that follows a negative binomial
  mat1 <- jitter(matrix(
    MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat1) <- paste0('CD', 1:ncol(mat1))

  mat2 <- jitter(matrix(
    MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat2) <- paste0('CD', 1:ncol(mat2))

  metadata <- data.frame(
    group = c('PB1', 'PB2'),
    row.names = c('mat1', 'mat2'),
    stringsAsFactors = FALSE)
}
