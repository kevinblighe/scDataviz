\name{downsampleByVar}

\alias{downsampleByVar}

\title{downsampleByVar}

\description{Downsample an input data-matrix based on variance.}

\usage{
  downsampleByVar(x, varianceFactor = 0.1)
}

\arguments{
  \item{x}{Input data-matrix. REQUIRED.}
  \item{varianceFactor}{Removes this proportion of variables based on
    lesser variance. DEFAULT = 0.1. OPTIONAL.}
}

\value{
A \code{\link{matrix}} object.
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
