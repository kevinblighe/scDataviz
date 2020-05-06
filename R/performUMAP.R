performUMAP <- function(
  indata,
  assay = 'scaled',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL,
  verbose = TRUE)
{
  if (is(indata, 'SingleCellExperiment')) {

    if (verbose) message('--input data class is SingleCellExperiment')

    if (!is.null(reducedDim)) {
        if (verbose) message('--input data is taken from \'', reducedDim,
          '\' dimensional reduction')
        if (verbose) message('--Dimensions to use: ', paste(dims, collapse = ', '))
        mat <- as.matrix(reducedDims(indata)[[reducedDim]][,dims])

        if (is.null(newDimName)) {
          newDimName <- paste0('UMAP_', reducedDim)
        }
    } else {
      if (verbose) message('--input data is taken from \'', assay, '\' assay slot')
      mat <- t(assay(indata, assay))

      if (is.null(newDimName)) {
        newDimName <- 'UMAP'
      }
    }

    if (verbose) message('--Performing UMAP...')
    if (is.null(useMarkers)) {
      u <- umap(mat)
    } else if (!is.null(useMarkers) && !is.null(reducedDim)) {
      warning('\'useMarkers\' and \'reducedDim\' are incompatible - ',
        'markers cannot be selected from a reduced dimensional ',
        'component in which they don\'t exist! Dimensions to use ',
        'for UMAP have already been chosen via the \'dims\' parameter')
      u <- umap(mat)
    } else {
      if (verbose) message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers])
    }

    colnames(u$layout) <- c('UMAP1', 'UMAP2')

    reducedDim(indata, newDimName) <- u$layout

    if (verbose) message('--Done')

    return(indata)

  } else {

    if (verbose) message('--input data class is ', class(indata))
    if (verbose) message('Note: all non-SingleCellExperiment objects will be ',
      'coerced to matrix')
    mat <- t(as.matrix(indata))

    if (verbose) message('--Performing UMAP...')
    if (is.null(useMarkers)) {
      u <- umap(mat)
    } else {
      if (verbose) message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers])
    }

    colnames(u$layout) <- c('UMAP1', 'UMAP2')

    if (verbose) message('--Done')

    return(u$layout)
  }
}
