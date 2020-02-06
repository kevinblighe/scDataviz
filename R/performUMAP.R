performUMAP <- function(
  sce,
  reducedDim = NULL,
  dims = c(1:20),
  newDimName = NULL,
  useMarkers = NULL)
{
  if (is.null(reducedDim)) {
    message('--input data is scaled expression levels')
    data <- t(assay(sce, 'scaled'))

    if (is.null(newDimName)) {
      newDimName <- 'UMAP'
    }
  } else if (reducedDim == 'PCA') {
    message('--input data is PC eigenvectors')
    data <- data.matrix(reducedDims(sce)[[reducedDim]][,dims])

    if (is.null(newDimName)) {
      newDimName <- 'UMAP_PCA'
    }
  }

  message('--Performing UMAP...')
  if (is.null(useMarkers)) {
    u <- umap(data)
  } else {
    u <- umap(data[,useMarkers])
  }

  colnames(u$layout) <- c('UMAP1', 'UMAP2')

  reducedDim(sce, newDimName) <- u$layout

  message('--Done')

  return(sce)
}
