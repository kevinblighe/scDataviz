#' Perform UMAP on an input data-frame or matrix, or \code{SingleCellExperiment} object, using the basic R implementation of UMAP.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'   data-frame or matrix, only the derived co-ordinates for the first two
#'   dimensions are returned. If a \code{SingleCellExperiment} object, UMAP is
#'   performed on the assay named by \code{assay}, and the co-ordinates for the
#'   first two dimensions are stored as a reduced dimension named by
#'   \code{reducedDim}.
#' @param config UMAP configuration settings
#' @param assay Name of the assay slot in \code{indata} from which data will be
#'   taken, assuming \code{indata} is a \code{SingleCellExperiment} object.
#' @param reducedDim A dimensional reduction / embedding stored within
#'   \code{indata}, e.g., PCA. If activated, UMAP will be performed on this object in
#'   place of the assay component specified by \code{assay}.
#' @param dims If 'reducedDim' is activated, the number of dimensions to use.
#' @param newDimName Name for the new dimensional embedding that will be produced.
#'   If nothing is selected for neither this nor \code{reducedDim}, then the new name
#'   will be \code{UMAP}. If nothing is selected for this, but, e.g., \code{PCA} is selected
#'   for \code{reducedDim}, then the new name will be \code{UMAP_PCA}.
#' @param useMarkers Before performing UMAP, subset the data for these markers.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Perform UMAP on an input data-frame or matrix, or \code{SingleCellExperiment} object, using the basic R implementation of UMAP.
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#'
#' performUMAP(mat)
#'
#' @import SingleCellExperiment
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom methods is
#'
#' @export
performUMAP <- function(
  indata,
  config = NULL,
  assay = 'scaled',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL,
  verbose = TRUE)
{
  if (is.null(config)) {
    config <- umap::umap.defaults
  }

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
      u <- umap(mat, config = config)
    } else if (!is.null(useMarkers) && !is.null(reducedDim)) {
      warning('\'useMarkers\' and \'reducedDim\' are incompatible - ',
        'markers cannot be selected from a reduced dimensional ',
        'component in which they don\'t exist! Dimensions to use ',
        'for UMAP have already been chosen via the \'dims\' parameter')
      u <- umap(mat, config = config)
    } else {
      if (verbose) message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers], config = config)
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
      u <- umap(mat, config = config)
    } else {
      if (verbose) message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers], config = config)
    }

    colnames(u$layout) <- c('UMAP1', 'UMAP2')

    if (verbose) message('--Done')

    return(u$layout)
  }
}
