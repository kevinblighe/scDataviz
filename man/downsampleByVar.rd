\name{downsampleByVar}

\alias{downsampleByVar}

\title{downsampleByVar}

\description{Downsample an input data-frame or matrix based on variance.}

\usage{
  downsampleByVar(x,
    varianceFactor = 0.1,
    verbose = TRUE)
}

\arguments{
  \item{x}{Input data-matrix. REQUIRED.}
  \item{varianceFactor}{Removes this proportion of variables based on
    lesser variance. DEFAULT = 0.1. OPTIONAL.}
  \item{verbose}{Boolean (TRUE / FALSE) to print messages to console or not.
    DEFAULT = TRUE. OPTIONAL.}
}

\value{
A \code{\link{matrix}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  # create random data that follows a negative binomial
  mat <- jitter(matrix(
    MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
    ncol = 20))

  downsampleByVar(mat, varianceFactor = 0.1)
}
