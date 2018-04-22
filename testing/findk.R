findk <- function(x, varianceFactor=5, FUNcluster, K.max=100, B=5)
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

	x <- x[apply(x, 1, FUN=function(x) sqrt(sum(x^2)))>EuclideanNormThreshold,]
	NoiseCorrected <- x
	NoiseCorrected[NoiseCorrected<BackgroundNoiseThreshold] <- 0
	x <- transFun(NoiseCorrected/asinhFactor)

	#Select out the top markers based on variance
	variances <- apply(x, 1, var)
	x <- x[order(variances, decreasing=TRUE),]
	x <- x[1:round((nrow(x)/varianceFactor),0),]

	#Compute the gap statistic for the data, using PAM clustering (builds clusters around a cell and distance from this)
	#K-means builds clusters around the mean of a group of cells
	GapStat <- clusGapKB(x, FUNcluster=FUNcluster, K.max=K.max, B=B)

	#Return the ideal number of clusters
	#	method can be one of "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax"
	idealclustersv1 <- maxSE(GapStat$Tab[,3], GapStat$Tab[,4], method="firstSEmax", SE.factor=1)
	idealclustersv2 <- maxSE(GapStat$Tab[,3], GapStat$Tab[,4], method="Tibs2001SEmax", SE.factor=1)
	idealclustersv3 <- maxSE(GapStat$Tab[,3], GapStat$Tab[,4], method="globalSEmax", SE.factor=1)
	idealclustersv4 <- maxSE(GapStat$Tab[,3], GapStat$Tab[,4], method="firstmax", SE.factor=1)
	idealclustersv5 <- maxSE(GapStat$Tab[,3], GapStat$Tab[,4], method="globalmax", SE.factor=1)

	cat(paste("\nFirst SE max=", idealclustersv1, sep=""))
	cat(paste("\nTibshirani et al. (2001)=", idealclustersv2, sep=""))
	cat(paste("\nDudoit & Fridlyand (2002)=", idealclustersv3, sep=""))
	cat(paste("\nFirst max=", idealclustersv4, sep=""))
	cat(paste("\nGlobal max=", idealclustersv5, sep=""))
	cat(paste("\n", sep=""))

	return(GapStat)
}
