markerEnrichment <- function(
  indata,
  meta = NULL,
  assay = 'scaled',
  metacluster,
  clusterAssign = metadata(indata)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  lowerPercentile = 5,
  upperPercentile = 5)
{
  iCellsPerCluster <- c()
  iTotalCells <- c()
  iPercentage <- c()
  NegativeMarkers <- c()
  PositiveMarkers <- c()

  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')
    metadata <- metadata(indata)
    data <- as.data.frame(t(as.matrix(assay(indata, assay))))

  } else {

    message('--input data class is ', class(indata))

    if (is.null(meta)) {
      stop('When the input data is a non-SingleCellExperiment object, ',
        '\'indata\' must relate to an expression matrix (cells as columns; ',
        'genes as rows), while \'meta\' must be non-NULL and relate to ',
        'metadata assocaited with this data.')
    } else if (!all(rownames(meta) == colnames(indata))) {
      stop('\'rownames(meta)\' must be equal to \'colnames(indata)\'')
    }

    metadata <- meta
    data <- as.data.frame(t(as.matrix(indata)))
  }

  data <- aggregate(data, list(clusterAssign), funcSummarise)
  data <- data[,-1]
  data <- apply(data, 2, scale, scale = FALSE)
  data <- t(data) / max(abs(range(data)))
  data <- rescale(data, c(-1,1))

  nclus <- length(unique(clusterAssign))

  res <- data.frame(row.names = seq(0, nclus-1))
  metares <- data.frame(row.names = seq(0, nclus-1))

  for (j in seq(0, nclus-1)) {
    iCellsPerCluster <- length(clusterAssign[clusterAssign == j])
    iTotalCells <- length(clusterAssign)
    iPercentage <- (iCellsPerCluster/iTotalCells) * 100

    NegativeMarkers <- names(which(data[,j] < (min(data[,j], na.rm = TRUE) + (((max(data[,j], na.rm = TRUE) - min(data[,j], na.rm = TRUE)) / 100) * lowerPercentile))))
    PositiveMarkers <- names(which(data[,j] > (max(data[,j], na.rm = TRUE) - (((max(data[,j], na.rm = TRUE) - min(data[,j], na.rm = TRUE)) / 100) * upperPercentile))))

    metaclusAbundance <- table(metadata[which(clusterAssign == j),metacluster])
    metaclusAbundance <- (metaclusAbundance / iCellsPerCluster) * 100
    studyvar <- table(metadata[which(clusterAssign == j),metacluster])

    res <- rbind(res,
      c(j,
        iCellsPerCluster,
        iTotalCells,
        iPercentage,
        paste(NegativeMarkers, '-', sep = '', collapse = ''),
        paste(PositiveMarkers, '+', sep = '', collapse = '')))
    colnames(res) <- c('Cluster', 'nCells','TotalCells','PercentCells','NegMarkers','PosMarkers')
    res$Cluster <- as.numeric(as.character(res$Cluster))
    res$nCells <- as.numeric(as.character(res$nCells))
    res$TotalCells <- as.numeric(as.character(res$TotalCells))
    res$PercentCells <- as.numeric(as.character(res$PercentCells))
    res$PosMarkers <- as.character(res$PosMarkers)
    res$NegMarkers <- as.character(res$NegMarkers)

    metares <- rbind(metares,
      c(metaclusAbundance, studyvar))
    colnames(metares) <- c(paste0('PerCent_', names(metaclusAbundance)), paste0('nCell_', names(studyvar)))
  }

  res$PosMarkers[res$PosMarkers == '+'] <- NA
  res$NegMarkers[res$NegMarkers == '-'] <- NA

  return(cbind(res, metares))
}
