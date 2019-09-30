processFCS <- function(
  files,
  metadata = NULL,
  bgNoiseThreshold = 1,
  euclideanNormThreshold = 1,
  transformation = FALSE,
  transFun = function (x) asinh(x),
  asinhFactor = 5,
  downsample = 10,
  colsDiscard = c('Time','Event_length','Center','Offset','Width','Residual','tSNE1','tSNE2'),
  colsRetain = NULL,
  newColnames = NULL)
{
  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as filelist
  if (!is.null(metadata)) {
    if(!identical(filelist, rownames(metadata))) {
      stop("'filelist' is not identical to 'rownames(metadata)'")
    }
  }

  # read in the data to a list
  samples <- list()
  samples <- lapply(filelist, function(x) exprs(read.FCS(x, transformation = FALSE)))
  names(samples) <- filelist

  # filter markers out
  if (!is.null(colsDiscard)) {
    samples <- lapply(
      samples,
      function(x) data.matrix(x[,-which(colnames(x) %in% colsDiscard)]))
  }

  # filter markers in
  if (!is.null(colsRetain)) {
    samples <- lapply(
      samples,
      function(x) data.matrix(x[,which(colnames(x) %in% colsRetain)]))
  }

  # rename markers
  if(!is.null(newColnames)) {
    samples <- lapply(
      samples,
      function(x) colnames(x) <- newColnames)
  }

  # transform
  if (!is.null(transformation)) {
    #Applying Euclidean norm
    samples <- lapply(
      samples,
      function(x) x[apply(x, 1, FUN = function(x) sqrt(sum(x^2))) > euclideanNormThreshold,])

    # noise correction
    for(i in 1:length(samples)) {
      x <- samples[[i]]
      x[x < bgNoiseThreshold] <- 0
      samples[[i]] <- x
    }

    # transform
    samples <- lapply(
      samples,
      function(x) transFun(x / asinhFactor))
  }

  #load function for downsampling
  if(!is.null(downsample)) {
    if (downsample > 0) {
      samples <- lapply(
        samples,
        function(x) downsampleByVar(x, varianceFactor = downsample))
    }
  }

  # is there metadata?
  names <- colnames(metadata)
  metanew <- list()
  if (!is.null(metadata)) {
    for (i in 1:length(samples)) {
      tmp <- data.frame(row.names = 1:nrow(samples[[i]]))
      for (j in 1:ncol(metadata)) {
        tmp <- cbind(tmp, rep(metadata[i,j], nrow(samples[[i]])))
      }
      metanew[[i]] <- tmp
    }

    metadata <- do.call(rbind, metanew)
    colnames(metadata) <- names
  }

  # combine all samples
  samples <- do.call(rbind, samples)
  rownames(samples) <- 1:nrow(samples)

  # these should be equal
  if (!is.null(metadata)) {
    if (nrow(metadata) != nrow(samples)) {
      stop(paste0('Metadata does not match expression data',
        ' - please check your input.'))
    }
  }

  # return a list object that contains expression data and
  # coresponding metadata
  return(list(
    expression = samples,
    metadata = metadata))
}