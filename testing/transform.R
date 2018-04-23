checkdistribution <- function (x, cex=1, main=main, breaks=breaks, ...)
{
	#Applying Euclidean norm
	x <- x[apply(x, 1, FUN=function(x) sqrt(sum(x^2)))>EuclideanNormThreshold,]

	NoiseCorrected <- x
	NoiseCorrected[NoiseCorrected<BackgroundNoiseThreshold] <- 0

	par(cex=cex)
	hist(transFun(NoiseCorrected/asinhFactor), main=main, breaks=breaks)
}
