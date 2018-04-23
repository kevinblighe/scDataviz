transform <- function(x, BackgroundNoiseThreshold, EuclideanNormThreshold, transFun, asinhFactor)
{
	#Applying Euclidean norm
	x <- x[apply(x, 1, FUN=function(x) sqrt(sum(x^2)))>EuclideanNormThreshold,]

	NoiseCorrected <- x
	NoiseCorrected[NoiseCorrected<BackgroundNoiseThreshold] <- 0

	x <- transFun(NoiseCorrected/asinhFactor)

	return(x)
}