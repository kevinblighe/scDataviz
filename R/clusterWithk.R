clusterWithk <- function(x, varianceFactor=5, FUNcluster, k=50, lowerPercentile=12.5, upperPercentile=12.5, p=0.05, ...)
{
	packageExists <- require(cluster)
	if(!packageExists)
	{
		stop( "Please install cluster first.", call.=FALSE)
	}

	packageExists <- require(fpc)
	if(!packageExists)
	{
		stop( "Please install fpc first.", call.=FALSE)
	}

	#Select out the top markers based on variance
	x <- downsampleByVar(x, varianceFactor)

	clus <- FUNcluster(x, k)

	clus <- clus[[1]]

	###

	iCellsPerCluster <- c()
	iTotalCells <- c()
	iPercentage <- c()
	NegativeMarkers <- c()
	PositiveMarkers <- c()

	#Count percentage of cells and determine which markers are expressed or not
	for (j in 1:k)
	{
		iCellsPerCluster <- length(clus$clustering[clus$clustering==j])
		iTotalCells <- length(clus$clustering)
		iPercentage <- (iCellsPerCluster/iTotalCells) * 100

		#Determine scaling factor
		DataMatrix <- apply(clus$medoids, 2, scale, scale=FALSE)
		DataMatrix <- t(DataMatrix) / max(abs(range(DataMatrix)))

		#Negative markers fall under 12.5%
		NegativeMarkers <- names(which(DataMatrix[,j] < (min(DataMatrix[,j]) + (((max(DataMatrix[,j]) - min(DataMatrix[,j])) / 100) * lowerPercentile))))

		#Positive markers rise above 50%
		PositiveMarkers <- names(which(DataMatrix[,j] > (max(DataMatrix[,j]) - (((max(DataMatrix[,j]) - min(DataMatrix[,j])) / 100) * upperPercentile))))

		print(paste(j, iPercentage, paste(PositiveMarkers, "+", sep="", collapse=""), paste(NegativeMarkers, "-", sep="", collapse=""), sep=", "))

		clus$iCellsPerCluster[j] <- iCellsPerCluster
		clus$iTotalCells[j] <- iTotalCells
		clus$iPercentage[j] <- iPercentage
		clus$NegativeMarkers[j] <- paste(NegativeMarkers, "-", sep="", collapse="")
		clus$PositiveMarkers[j] <- paste(PositiveMarkers, "+", sep="", collapse="")
	}

return(clus)
}
