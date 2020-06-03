#' Downsample an input data-frame or matrix based on variance.
#'
#' @param x Input data-matrix.
#' @param varianceFactor Removes this proportion of variables based on
#'   lesser variance.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Downsample an input data-frame or matrix based on variance.
#'
#' @return A \code{matrix} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
#'   ncol = 20))
#'
#' downsampleByVar(mat, varianceFactor = 0.1)
#'
#' @importFrom MASS rnegbin
#' @importFrom utils head
#' @importFrom matrixStats rowVars
#' 
#' @export
downsampleByVar <- function(x,
  varianceFactor = 0.1,
  verbose = TRUE)
{
  rowVars <- NULL

  vars <- rowVars(x)
  varorder <- order(vars, decreasing = TRUE)
  keep <- head(varorder, max(1, nrow(x)*(1-varianceFactor)))
  x <- x[keep,,drop=FALSE]

  return(x)
}
