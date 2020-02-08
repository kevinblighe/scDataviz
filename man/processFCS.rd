\name{processFCS}

\alias{processFCS}

\title{processFCS}

\description{Input, filter, normalise, and transform FCS expression data.}

\usage{
  processFCS(files,
  assayname = 'scaled',
  metadata = NULL,
  filter = TRUE,
  bgNoiseThreshold = 1,
  euclideanNormThreshold = 1,
  transformation = TRUE,
  transFun = function (x) asinh(x),
  asinhFactor = 5,
  downsample = 100000,
  downsampleVar = 0.1,
  colsDiscard = c('Time','Event_length','Center','Offset',
    'Width','Residual','tSNE1','tSNE2','BCKG'),
  colsRetain = NULL,
  newColnames = NULL)
}

\arguments{
  \item{files}{A vector of FCS files. REQUIRED.}
  \item{assay}{Name of the assay slot in sce from which data will be taken.
    DEFAULT = 'scaled'. OPTIONAL.}
  \item{metadata}{Metadata associated with the FCS files specified in
    'files'. A strict rule is enforced requiring that rownames(metadata)
    matches files in both name and order. DEFAULT = NULL. OPTIONAL.}
  \item{filter}{Boolean (TRUE / FALSE) to enable filtering (per sample)
    for background signal / noise. DEFAULT = TRUE. OPTIONAL.}
  \item{bgNoiseThreshold}{Threshold for background noise. Used when 
    'filter' == TRUE. DEFAULT = 1. OPTIONAL.}
  \item{euclideanNormThreshold}{Euclidean norm threshold for background
    noise. Used when 'filter' == TRUE. DEFAULT = 1. OPTIONAL.}
  \item{transformation}{Boolean (TRUE / FALSE) to enable data transformation
    after filtering. DEFAULT = TRUE. OPTIONAL.}
  \item{transFun}{The function to apply (per sample) for transformation. 
    Typically, for flow and mass cytometry, this is hyperbolic arc sine
    (asinh(x)). User can supply any function. DEFAULT = function (x) asinh(x).
    OPTIONAL.}
  \item{asinhFactor}{The factor to apply when transforming via asinh(). For
    flow cytometry, this is usually 150; for mass cytometry and CyTOF, it is
    5. Note that this is not used if the user has supplied their own function
    to 'transFun'. DEFAULT = 5. OPTIONAL.}
  \item{downsample}{Downsample to this number of random variables. This is
    perfromed on the final merged dataset, i.e., after all samples have been
      bound together. NULL to disable. ThiDEFAULT = 100000. OPTIONAL.}
  \item{downsampleVar}{Downsample based on variance. Removes this proportion of
    cells based on lesser variance. This is applied per sample. If user wishes
    to apply this globally on the final merged dataset, then set this to 0 and
    remove based on variance manually.  DEFAULT = 0.1. OPTIONAL.}
  \item{colsDiscard}{Columns to be removed from the final merged data. This
    names are literal and must match exactly. DEFAULT = c('Time','Event_length',
    'Center','Offset','Width','Residual','tSNE1','tSNE2','BCKG'). OPTIONAL.}
  \item{colsRetain}{Retain these columns only. This is the same as 'colsDiscard'
    but in reverse. Technically, it is possible to activate both 'colsDiscard'
    and 'colsRetain', but 'colsDiscard' will be executed first. DEFAULT = NULL.
    OPTIONAL.}
  \item{newColnames}{Assuming that you know the exact order of your final selected
    markers, rename these based on a vector passed as this argument. Please
    exercise caution when using this. DEFAULT = NULL. OPTIONAL.}
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
