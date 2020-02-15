\name{plotSignatures}

\alias{plotSignatures}

\title{plotSignatures}

\description{Find enriched markers per identified cluster and visualise these as a custom corrplot. plotSignatures first collapses your input data's expression profiles from the level of cells to the level of clusters based on a mathematical function specified by 'funcSummarise'. It then centers, independently per marker, the expression of each marker across all clusters. Finaly, data is scaled to be between -1 and +1.}

\usage{
  plotSignatures(
    indata,
    assay = 'scaled',
    clusterAssign = metadata(indata)[['Cluster']],
    funcSummarise = function(x) median(x, na.rm = TRUE),
    col = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100),
    labCex = 1.0,
    legendCex = 1.0,
    labDegree = 90)
}

\arguments{
  \item{indata}{A data-frame or matrix, or SingleCellExperiment object. If a
    data-frame or matrix, this should relate to expression data (cells as
    columns; genes as rows). If a SingleCellExperiment object, data will be
    extracted from an assay component named by 'assay'. REQUIRED.}
  \item{assay}{Name of the assay slot in 'indata' from which data will be
    taken, assuming 'indata' is a SingleCellExperiment object.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{clusterAssign}{A vector of cell-to-cluster assignments. This can be
    from any source but must align with your cells / variables. There is no
    check to ensure this when 'indata' is not a SingleCellExperiment object.
    DEFAULT = metadata(indata)[['Cluster']]. OPTIONAL.}
  \item{funcSummarise}{A mathematical function used to summarise expression
    per marker, per cluster. DEFAULT = function(x) median(x, na.rm = TRUE).
    OPTIONAL.}
  \item{col}{colorRampPalette to be used for shading low-to-high expression.
    DEFAULT = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100). OPTIONAL.}
  \item{labCex}{cex (size) of the main plot labels. DEFAULT = 1.0. OPTIONAL.}
  \item{legendCex}{cex (size) of the legend labels. DEFAULT = 1.0. OPTIONAL.}
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
  mat <- jitter(matrix(
    MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat) <- paste0('CD', 1:ncol(mat))
  rownames(mat) <- paste0('cell', 1:nrow(mat))

  u <- umap::umap(mat)$layout
  colnames(u) <- c('UMAP1','UMAP2')
  rownames(u) <- rownames(mat)
  clus <- clusKNN(u)

  plotSignatures(t(mat), clusterAssign = clus)
}
