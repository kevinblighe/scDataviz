plotSignatures <- function(clus=null, col=col, cexlab=1.0, cexlegend=1.0, labDegree=90)
{
	require(corrplot)

	DataMatrix <- clus$medoids
	DataMatrix <- DataMatrix[,order(colnames(DataMatrix))]

	#Center the medoids
	DataMatrix <- apply(DataMatrix, 2, scale, scale=FALSE)

	#Scale between -1 and +1
	DataMatrix <- t(DataMatrix) / max(abs(range(DataMatrix)))

	corrplot(t(DataMatrix), method="circle", order="original", addgrid.col="grey60", tl.col="black", col=col, cl.cex=cexlab, cl.pos="b", cl.ratio=0.4, cl.lim=c(-1,1), cl.length=3, tl.cex=cexlab, tl.srt=labDegree, mar=c(1,2,1,2))
	legend("top", bty="n", cex=cexlegend, title="", c("High", "Average", "Low"), fill=c(col[length(col)-10], col[length(col)/2], col[11]), horiz=TRUE)
}
