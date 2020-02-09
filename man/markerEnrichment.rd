\name{markerEnrichment}

\alias{markerEnrichment}

\title{markerEnrichment}

\description{Find enriched markers per identified cluster and calculate cluster abundances across these for metadata traits.}

\usage{
  markerEnrichment(sce,
  assay = 'scaled',
  metacluster,
  clusterVector = metadata(sce)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  lowerPercentile = 5,
  upperPercentile = 5)
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.}
  \item{assay}{Name of the assay slot in sce from which data will be taken.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{metacluster}{A column name from 'metadata(sce)' representing a trait
    over which metacluster abundances will be calculated. REQUIRED.}
  \item{clusterVector}{A vector of cell-to-cluster assignments. This can be
    from any source but ought to be taken from the metadata. DEFAULT =
    metadata(sce)[['Cluster']]. OPTIONAL.}
  \item{funcSummarise}{A mathematical function used to summarise expression
    per marker per cluster. DEFAULT = function(x) median(x, na.rm = TRUE).
    OPTIONAL.}
  \item{lowerPercentile}{Lower limit (%) of the expression range for selecting a
    signature marker for each cluster. DEFAULT = 5. OPTIONAL.}
  \item{upperPercentile}{Upper limit (%) of the expression range for selecting a
    signature marker for each cluster. DEFAULT = 5. OPTIONAL.}
}

\value{
A \code{\link{data.frame}} object.
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
