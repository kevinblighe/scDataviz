\name{performUMAP}

\alias{performUMAP}

\title{performUMAP}

\description{Perform UMAP on the scaled assay component of a SingleCellExperiment object, or a specified reduced dimensional component stored within this object.}

\usage{
  performUMAP(sce,
  reducedDim = NULL,
  dims = c(1:20),
  newDimName = NULL,
  useMarkers = NULL)
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.}
  \item{reducedDim}{A reduced dimensional component stored within 'sce',
    e.g., PCA. If activated, UMAP will be performed on this object.
    and not the scaled assay component of 'sce'. DEFAULT = NULL. OPTIONAL.}
  \item{dims}{If 'reducedDim' isa activated, the number of dimensions to
    use. DEFAULT = c(1:20). OPTIONAL.}
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
  # not run
}
