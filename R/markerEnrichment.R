markerEnrichment <- function(
  data,
  lowerPercentile = 5,
  upperPercentile = 5)
{ 
  iCellsPerCluster <- c()
  iTotalCells <- c()
  iPercentage <- c()
  NegativeMarkers <- c()
  PositiveMarkers <- c()
  res <- data.frame()

  df <- data$expression
  df <- aggregate(data.matrix(data$expression), data$nnc, mean)
  df <- df[,-1]
  df <- apply(df, 2, scale, scale = FALSE)
  df <- t(df) / max(abs(range(df)))
  df <- scales::rescale(df, c(-1,1))

  clus <- data$nnc[,1]

  #Count percentage of cells and determine which markers are expressed or not
  nclus <- length(unique(clus))
  for (j in 0:(nclus-1)) {
    iCellsPerCluster <- length(clus[clus == j])
    iTotalCells <- length(clus)
    iPercentage <- (iCellsPerCluster/iTotalCells) * 100
    NegativeMarkers <- names(which(df[,j] < (min(df[,j]) + (((max(df[,j]) - min(df[,j])) / 100) * lowerPercentile))))
    PositiveMarkers <- names(which(df[,j] > (max(df[,j]) - (((max(df[,j]) - min(df[,j])) / 100) * upperPercentile))))
    res <- rbind(
      res,
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
  }

  return(res)
}
