\name{markerEnrichment}

\alias{markerEnrichment}

\title{markerEnrichment}

\description{Find enriched markers per identified cluster and calculate cluster abundances across these for metadata traits.}

\usage{
  markerEnrichment(sce,
  metacluster,
  clusterVector = metadata(sce)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  lowerPercentile = 5,
  upperPercentile = 5)
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.}
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
  # not run
}
