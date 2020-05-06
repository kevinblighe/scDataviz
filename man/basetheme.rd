\name{basetheme}

\alias{basetheme}

\title{basetheme}

\description{Package-wide, non-user function used to set a base ggplot2 theme.}

\usage{
  basetheme(
    titleLabSize,
    subtitleLabSize,
    captionLabSize,
    axisLabSize,
    xlabAngle,
    xlabhjust,
    xlabvjust,
    ylabAngle,
    ylabhjust,
    ylabvjust,
    legendPosition,
    legendLabSize)
  }

\arguments{
  \item{titleLabSize}{Size of plot title. REQUIRED.}
  \item{subtitleLabSize}{Size of plot subtitle. REQUIRED.}
  \item{captionLabSize}{Size of plot caption. REQUIRED.}
  \item{axisLabSize}{Size of x- and y-axis labels. REQUIRED.}
  \item{xlabAngle}{Rotation angle of x-axis labels. REQUIRED.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. REQUIRED.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. REQUIRED.}
  \item{ylabAngle}{Rotation angle of y-axis labels. REQUIRED.}
  \item{ylabhjust}{Horizontal adjustment of y-axis labels. REQUIRED.}
  \item{ylabvjust}{Vertical adjustment of y-axis labels. REQUIRED.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
    'none'). REQUIRED.}
  \item{legendLabSize}{Size of plot legend text. REQUIRED.}
}

\value{
A \code{\link{list}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  # create a theme
  th <- basetheme(
    titleLabSize = 16,
    subtitleLabSize = 12,
    captionLabSize = 12,
    axisLabSize = 16,
    xlabAngle = 0,
    xlabhjust = 0.5,
    xlabvjust = 0.5,
    ylabAngle = 0,
    ylabhjust = 0.5,
    ylabvjust = 0.5,
    legendPosition = 'none',
    legendLabSize = 12)
}
