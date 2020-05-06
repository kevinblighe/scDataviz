\name{markerEnrichment}

\alias{markerEnrichment}

\title{markerEnrichment}

\description{Find enriched markers per identified cluster and calculate cluster abundances across these for metadata traits. markerEnrichment first collapses your input data's expression profiles from the level of cells to the level of clusters based on a mathematical function specified by 'funcSummarise'. It then centers, independently per marker, the expression of each marker across all clusters. Finaly, data is scaled to be between -1 and +1.}

\usage{
  markerEnrichment(
    indata,
    meta = NULL,
    assay = 'scaled',
    metacluster,
    clusterAssign = metadata(indata)[['Cluster']],
    funcSummarise = function(x) median(x, na.rm = TRUE),
    lowerPercentile = 5,
    upperPercentile = 5,
    verbose = TRUE)
}

\arguments{
  \item{indata}{A data-frame or matrix, or SingleCellExperiment object. If a
    data-frame or matrix, this should relate to expression data (cells as
    columns; genes as rows). If a SingleCellExperiment object, data will be
    extracted from an assay component named by 'assay'. REQUIRED.}
  \item{meta}{If 'indata' is a non-SingleCellExperiment object, 'meta' must be
    activated and relate to a data-frame of metadata that aligns with the columns
    of 'indata', and that also contains a column name specified by 'metacluster'.
    DEFAULT = NULL. OPTIONAL.}
  \item{assay}{Name of the assay slot in 'indata' from which data will be
    taken, assuming 'indata' is a SingleCellExperiment object.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{metacluster}{A column name from the provided metadata representing a
    trait over which metacluster abundances will be calculated. REQUIRED.}
  \item{clusterAssign}{A vector of cell-to-cluster assignments. This can be
    from any source but must align with your cells / variables. There is no
    check to ensure this when 'indata' is not a SingleCellExperiment object.
    DEFAULT = metadata(indata)[['Cluster']]. OPTIONAL.}
  \item{funcSummarise}{A mathematical function used to summarise expression
    per marker per cluster. DEFAULT = function(x) median(x, na.rm = TRUE).
    OPTIONAL.}
  \item{lowerPercentile}{Lower limit (%) of the scaled expression range for
    selecting a signature marker for each cluster. DEFAULT = 5. OPTIONAL.}
  \item{upperPercentile}{Upper limit (%) of the scaled expression range for
    selecting a signature marker for each cluster. DEFAULT = 5. OPTIONAL.}
  \item{verbose}{Boolean (TRUE / FALSE) to print messages to console or not.
    DEFAULT = TRUE. OPTIONAL.}
}

\value{
A \code{\link{data.frame}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  # create random data that follows a negative binomial
  mat <- jitter(matrix(
    MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat) <- paste0('CD', 1:ncol(mat))
  rownames(mat) <- paste0('cell', 1:nrow(mat))

  u <- umap::umap(mat)$layout
  colnames(u) <- c('UMAP1','UMAP2')
  rownames(u) <- rownames(mat)
  clus <- clusKNN(u)

  metadata <- data.frame(
    group = c(rep('PB1', 25), rep('PB2', 25)),
    row.names = rownames(u))

  markerEnrichment(t(mat), meta = metadata,
    metacluster = 'group', clusterAssign = clus)
}
