facsplot <- function(marker1, marker2, df, bins=200, main="FACS plot", xlim=xlim, ylim=ylim, x1, x2, y1, y2, cex=1, colramp=rf, ...)
{
	packageExists <- require(hexbin)
	if(!packageExists)
	{
		stop( "Please install hexbin first.", call.=FALSE)
	}

	packageExists <- require(lattice)
	if(!packageExists)
	{
		stop( "Please install lattice first.", call.=FALSE)
	}

	x <- df[,marker1]
	y <- df[,marker2]

	x[x==0] <- NA
	x <- log2(x)

	y[y==0] <- NA
	y <- log2(y)

	par=mar=c(5,5,5,5, cex=0.8)

	hexbinplot(y ~ x, data=df, aspect="1", xbins=bins, xlab=marker1, ylab=marker2, xlim=xlim, ylim=ylim, cex.labels=1.0, cex.title=1.0, colramp=rf,
    		panel=function(x, y, ...)
    		{
        		panel.hexbinplot(x, y, ...)
        		panel.abline(v=c(x1,x2), h=c(y1,y2), col="black", lwd=4, lty=5)
    		})
}
