performUMAP <- function(
  data,
  useMarkers = NULL)
{
  u <- umap(data$expression)

  data$layout <- u$layout
  colnames(data$layout) <- c('UMAP1', 'UMAP2')

  return(data)
}