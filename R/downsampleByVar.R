downsampleByVar <- function(x,
  varianceFactor = 0.1,
  verbose = TRUE)
{
  rowVars <- NULL

  vars <- rowVars(x)
  if (verbose) message('--removing the lower ', varianceFactor * 100,
      '% of cells based on variance')
  varorder <- order(vars, decreasing = TRUE)
  keep <- head(varorder, max(1, nrow(x)*(1-varianceFactor)))
  x <- x[keep,,drop=FALSE]

  return(x)
}
