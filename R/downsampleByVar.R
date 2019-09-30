downsampleByVar <- function(x, varianceFactor = 5)
{
  variances <- apply(x, 1, var)
  x <- x[order(variances, decreasing=TRUE),]
  x <- x[1:round((nrow(x)/ceiling(varianceFactor)), 0),]
  
  return(x)
}