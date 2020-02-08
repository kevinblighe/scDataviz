\name{performUMAP}

\alias{performUMAP}

\title{performUMAP}

\description{Perform UMAP on the scaled assay component of a SingleCellExperiment object, or a specified reduced dimensional component stored within this object.}

\usage{
  performUMAP(sce,
  assay = 'scaled',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL)
}

\arguments{
  \item{sce}{A SingleCellExperiment object. REQUIRED.}
  \item{assay}{Name of the assay slot in sce from which data will be taken.
    DEFAULT = 'scaled'. OPTIONAL.}
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
  # not run
}
