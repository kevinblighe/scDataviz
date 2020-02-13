\name{performUMAP}

\alias{performUMAP}

\title{performUMAP}

\description{Perform UMAP using the basic R implementation of this.}

\usage{
  performUMAP(
    indata,
    assay = 'scaled',
    reducedDim = NULL,
    dims = seq_len(20),
    newDimName = NULL,
    useMarkers = NULL)
}

\arguments{
  \item{indata}{A data-frame or matrix, or SingleCellExperiment object. If a
    data-frame or matrix, only the derived co-ordinates for the first two
    dimensions are returned. If a SingleCellExperiment object, UMAP is
    performed on the assay named by 'assay', and the co-ordinates for the
    first two dimensions are stored as a reduced dimension named by
    'reducedDim'. REQUIRED.}
  \item{assay}{Name of the assay slot in 'indata' from which data will be
    taken, assuming 'indata' is a SingleCellExperiment object. DEFAULT = 'scaled'.
    OPTIONAL.}
  \item{reducedDim}{A reduced dimensional component stored within 'sce',
    e.g., PCA. If activated, UMAP will be performed on this object.
    and not the scaled assay component of 'sce'. DEFAULT = NULL. OPTIONAL.}
  \item{dims}{If 'reducedDim' isa activated, the number of dimensions to
    use. DEFAULT = seq_len(20). OPTIONAL.}
  \item{newDimName}{Name for the new dimensional component. If nothing is
    selected for this or 'reducedDim', then the new name will be 'UMAP'. If
    nothing is selected for this but 'PCA' is selected as 'reducedDim', then
    the new name will be 'UMAP_PCA'. DEFAULT = NULL. OPTIONAL.}
  \item{useMarkers}{Before performing UMAPsubset the data for these markers.
    DEFAULT = NULL. OPTIONAL.}
}

\value{
A \code{\link{SingleCellExperiment}} object.
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

  performUMAP(mat)
}
