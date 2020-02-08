markerEnrichment <- function(
  sce,
  assay = 'scaled',
  metacluster,
  clusterVector = metadata(sce)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  lowerPercentile = 5,
  upperPercentile = 5)
{
  iCellsPerCluster <- c()
  iTotalCells <- c()
  iPercentage <- c()
  NegativeMarkers <- c()
  PositiveMarkers <- c()

  metadata = metadata(sce)

  data <- as.data.frame(t(as.matrix(assay(sce, assay))))
  data <- aggregate(data, list(clusterVector), funcSummarise)
  data <- data[,-1]
  data <- apply(data, 2, scale, scale = FALSE)
  data <- t(data) / max(abs(range(data)))
  data <- rescale(data, c(-1,1))

  nclus <- length(unique(clusterVector))

  res <- data.frame(row.names = 0:(nclus-1))
  metares <- data.frame(row.names = 0:(nclus-1))

  for (j in 0:(nclus-1)) {
    iCellsPerCluster <- length(clusterVector[clusterVector == j])
    iTotalCells <- length(clusterVector)
    iPercentage <- (iCellsPerCluster/iTotalCells) * 100

    NegativeMarkers <- names(which(data[,j] < (min(data[,j]) + (((max(data[,j]) - min(data[,j])) / 100) * lowerPercentile))))
    PositiveMarkers <- names(which(data[,j] > (max(data[,j]) - (((max(data[,j]) - min(data[,j])) / 100) * upperPercentile))))

    metaclusAbundance <- table(metadata[which(clusterVector == j),metacluster])
    metaclusAbundance <- (metaclusAbundance / iCellsPerCluster) * 100
    studyvar <- table(metadata[which(clusterVector == j),metacluster])

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
