performUMAP <- function(
  sce,
  assay = 'scaled',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL)
{
  if (!is.null(reducedDim)) {
      message('--input data is taken from \'', reducedDim, '\' dimensional reduction')
      message('--Dimensions to use: ', paste(dims, collapse = ', '))
      data <- data.matrix(reducedDims(sce)[[reducedDim]][,dims])

      if (is.null(newDimName)) {
        newDimName <- paste0('UMAP_', reducedDim)
      }
  } else {
    message('--input data is taken from \'', assay, '\' assay slot')
    data <- t(assay(sce, assay))

    if (is.null(newDimName)) {
      newDimName <- 'UMAP'
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
