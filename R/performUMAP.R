performUMAP <- function(
  sce,
  useMarkers = NULL,
  stratify = NULL)
{
  data <- t(assay(sce, 'scaled'))

  if (is.null(stratify)) {
    message('--Performing UMAP for all data combined')
    if (is.null(useMarkers)) {
      u <- umap(data)
    } else {
      u <- umap(data[,useMarkers])
    }

    colnames(u$layout) <- c('UMAP1', 'UMAP2')

    reducedDims(sce) <- list(UMAP = u$layout)

    return(sce)
  } else {
    groups <- unique(metadata(sce)[,stratify])

    dims <- list()

    for (i in 1:length(groups)) {
      message(paste0('--Performing UMAP for ', groups[i], ' (', i, '/', length(groups), ')'))
      keep <- which(metadata(sce)[,stratify] == groups[i])

      if (is.null(useMarkers)) {
        u <- umap(data[keep,])
      } else {
        u <- umap(data[keep,useMarkers])
      }

      colnames(u$layout) <- c('UMAP1', 'UMAP2')

      dims[[i]] <- u$layout
      names(dims)[i] <- groups[i]
    }

    reducedDims(sce) <- dims
    return(sce)
  }
}
