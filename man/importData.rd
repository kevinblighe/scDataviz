\name{importData}

\alias{importData}

\title{importData}

\description{Import a data-frame or matrix, and associated metadata, to a SingleCellExperiment object.}

\usage{
  importData(
    mat,
    assayname,
    metadata = NULL,
    downsampleVar = NULL,
    verbose = TRUE)
}

\arguments{
  \item{mat}{A data-frame or matrix of expression values. Data-frames will be
    coerced to matrices. REQUIRED.}
  \item{assayname}{Name of the SingleCellExperiment assay slot in which
    data will be stored. REQUIRED.}
  \item{metadata}{Metadata associated with the data contained in 'mat'. A
    strict rule is enforced requiring that rownames(metadata) ==
    rownames(mat). DEFAULT = NULL. OPTIONAL.}
  \item{downsampleVar}{Downsample based on variance. Removes this proportion of
    variables (rows) based on lesser variance. This is applied on a per sample
    basis. If user wishes to apply this globally on the final merged dataset,
    then set this to 0 and remove based on variance manually after object
    creation. DEFAULT = NULL. OPTIONAL.}
  \item{verbose}{Boolean (TRUE / FALSE) to print messages to console or not.
    DEFAULT = TRUE. OPTIONAL.}
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
    MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat) <- paste0('CD', 1:ncol(mat))
  rownames(mat) <- paste0('cell', 1:nrow(mat))

  metadata <- data.frame(
    group = rep('A', nrow(mat)),
    row.names = rownames(mat),
    stringsAsFactors = FALSE)

  sce <- importData(mat,
    assayname = 'normcounts',
    metadata = metadata)
}
