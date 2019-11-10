processData <- function(
  mat,
  metadata = NULL,
  filter = TRUE,
  bgNoiseThreshold = 1,
  euclideanNormThreshold = 1,
  transformation = TRUE,
  transFun = function (x) asinh(x),
  asinhFactor = 5,
  downsample = 0.1,
  colsDiscard = c('Time','Event_length','Center','Offset','Width','Residual','tSNE1','tSNE2','BCKG'),
  colsRetain = NULL,
  newColnames = NULL)
{
  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as input list names
  if (!is.null(metadata)) {
    if(!identical(names(mat), rownames(metadata))) {
      stop("'mat' list name order is not identical to 'rownames(metadata)'")
    }
  }

  # carry over to same workflow as 'processFCS'
  samples <- mat

  # filter markers out
  if (!is.null(colsDiscard)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsDiscard)) > 0) {
        x[,-which(colnames(x) %in% colsDiscard)]} else {return(x)})
  }

  # filter markers in
  if (!is.null(colsRetain)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsDiscard)) > 0) {
        x[,which(colnames(x) %in% colsDiscard)]} else {return(x)})
  }

  # rename markers
  if(!is.null(newColnames)) {
    for(i in 1:length(samples)) {
      colnames(samples[[i]]) <- newColnames
    }
  }

  # filter
  if (filter == TRUE) {
    message('--filtering background / noise')

    # Euclidean norm
    samples <- lapply(
      samples,
      function(x) x[apply(x, 1, FUN = function(x) sqrt(sum(x^2))) > euclideanNormThreshold,])

    # noise correction
    for(i in 1:length(samples)) {
      x <- samples[[i]]
      x[x < bgNoiseThreshold] <- 0
      samples[[i]] <- x
    }
  }

  # transform
  if (transformation == TRUE) {
    message('--transforming data')
    samples <- lapply(
      samples,
      function(x) transFun(x / asinhFactor))
  }

  # load function for downsampling
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

  # return a SingleCellExperiment object
  ret <- SingleCellExperiment(
    assays = list(scaled = samples))
  rowData(ret) <- metadata
  return(ret)
}
