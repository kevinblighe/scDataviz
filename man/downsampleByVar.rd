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
 # not run
}
