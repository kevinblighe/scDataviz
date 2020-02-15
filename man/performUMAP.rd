\name{performUMAP}

\alias{performUMAP}

\title{performUMAP}

\description{Perform UMAP on an input data-frame or matrix, or SingleCellExperiment object, using the basic R implementation of UMAP.}

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
    taken, assuming 'indata' is a SingleCellExperiment object. DEFAULT =
    'scaled'. OPTIONAL.}
  \item{reducedDim}{A dimensional reduction / embedding stored within
    'indata', e.g., PCA. If activated, UMAP will be performed on this object in
    place of the assay component specified by 'assay'. DEFAULT = NULL. OPTIONAL.}
  \item{dims}{If 'reducedDim' is activated, the number of dimensions to use.
    DEFAULT = seq_len(20). OPTIONAL.}
  \item{newDimName}{Name for the new dimensional embedding that will be produced.
    If nothing is selected for neither this nor 'reducedDim', then the new name
    will be 'UMAP'. If nothing is selected for this, but, e.g., 'PCA' is selected
    for 'reducedDim', then the new name will be 'UMAP_PCA'. DEFAULT = NULL.
    OPTIONAL.}
  \item{useMarkers}{Before performing UMAP, subset the data for these markers.
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
