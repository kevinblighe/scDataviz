#' Input, filter, normalise, and transform FCS expression data.
#'
#' @param files A vector of FCS files.
#' @param assayname Name of the assay slot in which data will be stored.
#' @param metadata Metadata associated with the FCS files specified in
#'   'files'. A strict rule is enforced requiring that \code{rownames(metadata)}
#'   matches files in both name and order.
#' @param filter Boolean (TRUE / FALSE) to enable filtering (per sample)
#'   for background signal / noise.
#' @param bgNoiseThreshold Threshold for background noise. Used when 
#'   \code{filter == TRUE}.
#' @param euclideanNormThreshold Euclidean norm threshold for background
#'   noise. Used when \code{filter == TRUE}.
#' @param transformation Boolean (TRUE / FALSE) to enable data transformation
#'   after filtering.
#' @param transFun The function to apply (per sample) for transformation. 
#'   Typically, for flow and mass cytometry, this is hyperbolic arc sine
#'   (\code{asinh(x)}). User can supply any function.
#' @param asinhFactor The factor to apply when transforming via \code{asinh()}. For
#'   flow cytometry, this is usually 150; for mass cytometry and CyTOF, it is
#'   5. Note that this is not used if the user has supplied their own function
#'   to \code{transFun}.
#' @param downsample Downsample to this number of random variables. This is
#'   perfromed on the final merged dataset, i.e., after all samples have been
#'   bound together. NULL to disable.
#' @param downsampleVar Downsample based on variance. Removes this proportion of
#'   cells based on lesser variance. This is applied per sample. If user wishes
#'   to apply this globally on the final merged dataset, then set this to 0 and
#'   remove based on variance manually.
#' @param colsDiscard Columns to be removed from the final merged data. These
#'   names are literal and must match exactly.
#' @param colsRetain Retain these columns only. This is the same as \code{colsDiscard}
#'   but in reverse. Technically, it is possible to activate both \code{colsDiscard}
#'   and \code{colsRetain}, but \code{colsDiscard} will be executed first.
#' @param newColnames Assuming that you know the exact order of your final selected
#'   markers, rename these based on a vector passed as this argument. Please
#'   exercise caution when using this.
#' @param emptyValue boolean (taken from ?flowCore::read.FCS indicating whether or
#'   not we allow an empty value for keyword values in TEXT segment.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Input, filter, normalise, and transform FCS expression data.
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat1 <- jitter(matrix(
#'   MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat1) <- paste0('CD', 1:ncol(mat1))
#' rownames(mat1) <- paste0('cell', 1:nrow(mat1))
#'
#' mat2 <- jitter(matrix(
#'   MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat2) <- paste0('CD', 1:ncol(mat2))
#' rownames(mat2) <- paste0('cell', 1:nrow(mat2))
#'
#' metadata <- data.frame(
#'   group = c('PB1', 'PB2'),
#'   row.names = c('mat1', 'mat2'),
#'   stringsAsFactors = FALSE)
#'
#' @import SingleCellExperiment
#' 
#' @importFrom MASS rnegbin
#' @importFrom flowCore read.FCS exprs
#' @importFrom S4Vectors metadata<-
#'
#' @export
processFCS <- function(
  files,
  assayname = 'scaled',
  metadata = NULL,
  filter = TRUE,
  bgNoiseThreshold = 1,
  euclideanNormThreshold = 1,
  transformation = TRUE,
  transFun = function (x) asinh(x),
  asinhFactor = 5,
  downsample = 100000,
  downsampleVar = 0.1,
  colsDiscard = c('Time','Event_length','Center','Offset','Width',
    'Residual','tSNE1','tSNE2','BCKG'),
  colsRetain = NULL,
  newColnames = NULL,
  emptyValue = TRUE,
  verbose = TRUE)
{

  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as filelist
  if (!is.null(metadata)) {
    if(!identical(files, rownames(metadata))) {
      stop("'filelist' is not identical to 'rownames(metadata)'")
    }
  }

  # read in the data to a list
  samples <- list()
  samples <- lapply(files,
    function(x) exprs(read.FCS(x, transformation = FALSE, emptyValue = emptyValue)))
  names(samples) <- files

  # filter markers out
  if (!is.null(colsDiscard)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsDiscard))) {
        x[,-which(colnames(x) %in% colsDiscard)]} else {return(x)})
  }

  # filter markers in
  if (!is.null(colsRetain)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsRetain))) {
        x[,which(colnames(x) %in% colsRetain)]} else {return(x)})
  }

  # rename markers
  if(!is.null(newColnames)) {
    for(i in seq_len(length(samples))) {
      colnames(samples[[i]]) <- newColnames
    }
  }

  # filter
  if (filter) {
    if (verbose) message('--filtering background / noise')

    # Euclidean norm
    samples <- lapply(
      samples,
      function(x)
        x[apply(x, 1, FUN = function(x) sqrt(sum(x^2))) > euclideanNormThreshold,])

    # noise correction
    for(i in seq_len(length(samples))) {
      x <- samples[[i]]
      x[x < bgNoiseThreshold] <- 0
      samples[[i]] <- x
    }
  }

  # transform
  if (transformation) {
    if (verbose) message('--transforming data')
    samples <- lapply(
      samples,
      function(x) transFun(x / asinhFactor))
  }

  # load function for downsampling based on variance
  if(!is.null(downsampleVar)) {
    if (downsampleVar > 0) {
      if (verbose) message('--removing the lower ', downsampleVar * 100,
        '% of cells based on variance')

      samples <- lapply(
        samples,
        function(x) downsampleByVar(x, varianceFactor = downsampleVar,
          verbose = verbose))
    }
  }

  # is there metadata?
  names <- colnames(metadata)
  metanew <- list()
  if (!is.null(metadata)) {
    for (i in seq(length(samples))) {
      tmp <- data.frame(row.names = seq_len(nrow(samples[[i]])))
      for (j in seq_len(ncol(metadata))) {
        tmp <- cbind(tmp, rep(metadata[i,j], nrow(samples[[i]])))
      }
      metanew[[i]] <- tmp
    }

    metadata <- do.call(rbind, metanew)
    colnames(metadata) <- names
    rownames(metadata) <- paste0('cell', seq_len(nrow(metadata)))
  }

  # combine all samples
  samples <- do.call(rbind, samples)
  rownames(samples) <- paste0('cell', seq_len(nrow(samples)))

  # downsample
  if (!is.null(downsample)) {
    if (downsample > nrow(samples)) {
      warning('Cannot downsample to ', downsample, ' number of variables as',
        ' there are ', nrow(samples), ' variables currently in the merged ',
        'dataset.')
      if (verbose) message('--skipping downsampling')
    } else {
      if (verbose) message('--downsampling to ', downsample, ' variables.')
      idx <- sample(seq(nrow(samples)), downsample)
      samples <- samples[idx,]
      metadata <- metadata[idx,]

      rownames(metadata) <- paste0('cell', seq_len(nrow(metadata)))
      rownames(samples) <- paste0('cell', seq_len(nrow(samples)))
    }
  }

  # these should be equal
  if (!is.null(metadata)) {
    if (nrow(metadata) != nrow(samples)) {
      stop(paste0('Metadata does not match expression data',
        ' - please check your input.'))
    }
  }

  # return a SingleCellExperiment object
  ret <- list(t(samples))
  names(ret)[1] <- assayname
  ret <- SingleCellExperiment(
    assays = ret)
  metadata(ret) <- metadata
  return(ret)
}
