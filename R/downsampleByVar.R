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
