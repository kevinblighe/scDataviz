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
