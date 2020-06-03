#' Import a data-frame or matrix, and associated metadata, to a \code{SingleCellExperiment} object.
#'
#' @param mat A data-frame or matrix of expression values. Data-frames will be
#'   coerced to matrices.
#' @param assayname Name of the \code{SingleCellExperiment} assay slot in which
#'   data will be stored.
#' @param metadata Metadata associated with the data contained in 'mat'. A
#'   strict rule is enforced requiring that \code{rownames(metadata) ==
#'   rownames(mat)}.
#' @param downsampleVar Downsample based on variance. Removes this proportion of
#'   variables (rows) based on lesser variance. This is applied on a per sample
#'   basis. If user wishes to apply this globally on the final merged dataset,
#'   then set this to 0 and remove based on variance manually after object
#'   creation.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Import a data-frame or matrix, and associated metadata, to a \code{SingleCellExperiment} object.
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#' rownames(mat) <- paste0('cell', 1:nrow(mat))
#'
#' metadata <- data.frame(
#'   group = rep('A', nrow(mat)),
#'   row.names = rownames(mat),
#'   stringsAsFactors = FALSE)
#'
#' sce <- importData(mat,
#'   assayname = 'normcounts',
#'   metadata = metadata)
#'
#' @import SingleCellExperiment
#'
#' @importFrom MASS rnegbin
#' @importFrom S4Vectors metadata<-
#' 
#' @export
importData <- function(
  mat,
  assayname,
  metadata = NULL,
  downsampleVar = NULL,
  verbose = TRUE)
{

  # avoid attempting to coerce S4 matrices into full matrices.
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }

  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as rownames(mat)
  if (!is.null(metadata)) {
    if(!identical(rownames(mat), rownames(metadata))) {
      stop("'rownames(mat)' is not identical to 'rownames(metadata)'")
    }
  }

  # remove lower portion of variables based on variation
  vars <- rowVars(mat)
  if (!is.null(downsampleVar)) {
    if (verbose) message('-- removing the lower ', downsampleVar * 100,
      '% of variables based on variance')

    varorder <- order(vars, decreasing = TRUE)

    keep <- head(varorder, max(1, nrow(mat)*(1-downsampleVar)))

    mat <- mat[keep,,drop=FALSE]

    if (!is.null(metadata)) {
      metadata <- metadata[keep,,drop = FALSE]
    }

    vars <- vars[keep]
  }

  # if no metadata supplied, then create an empty data-frame
  # with just rownames
  if (is.null(metadata)) {
    metadata <- data.frame(row.names = rownames(mat))
  }

  # these should be equal at this point
  if (nrow(metadata) != nrow(mat)) {
    stop(paste0('Metadata does not match expression data',
      ' - please check your input.'))
  }

  # return a SingleCellExperiment object
  ret <- list(t(mat))
  names(ret)[1] <- assayname
  ret <- SingleCellExperiment(
    assays = ret)
  metadata(ret) <- metadata
  return(ret)
}
