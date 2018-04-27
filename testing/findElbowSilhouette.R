findElbowSilhouette <- function(x, varianceFactor=5, K.max=100, B=5)
{
	packageExists <- require(factoextra)
	if(!packageExists)
	{
		stop( "Please install factoextra first.", call.=FALSE)
	}

	packageExists <- require(NbClust)
	if(!packageExists)
	{
		stop( "Please install NbClust first.", call.=FALSE)
	}

	#Select out the top markers based on variance
	x <- downsampleByVar(x, varianceFactor)

	elbow <- fviz_nbclust(x, FUNcluster=CustomPAM, method="wss", diss=NULL, k.max=K.max, nboot=B)

	sil <- fviz_nbclust(x, FUNcluster=CustomPAM, method="silhouette", diss=NULL, k.max=K.max, nboot=B, print.summary=FALSE)

	return(list(elbow, sil))
}