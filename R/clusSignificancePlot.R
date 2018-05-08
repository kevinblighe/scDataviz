#Function:	clusSignificancePlot
#Requires:	lattice, latticeExtra, RColorBrewer
clusSignificancePlot <- function(df, labcex, strPalette, iNumColours, boolReverse, strTitle)
{
  #Required packages
  require(lattice)
  require(latticeExtra)
  require(RColorBrewer)

  #Check to see if everything is numeric; if not, print a warning message
  for (i in 1:ncol(df))
  {
    if(!is.numeric(df[,i]))
    {
      print(paste("Warning: ", x[i], " is not numeric - please check the source data as everything will be converted to a matrix", sep=""))
    }
  }

  #Save the p-values
  pvals <- round(data.matrix(df),2)

  #Convert to negative log (base 10)
  df <- -log10(data.matrix(df))
  df[df==Inf] <- 1

  #Determine max and min values in order to define the range
  max <- max(df)
  min <- min(df)
  iUpperRange <- max+1
  iLowerRange <- 0

  #Define the colour scheme/palette
  if (boolReverse==TRUE)
  {
    cols <- colorRampPalette(rev(brewer.pal(iNumColours, strPalette)))
  }
  else
  {
    cols <- colorRampPalette(brewer.pal(iNumColours, strPalette))
  }

  #Define a panel function for adding labels
  #Labels are passed with z as a third dimension
  labels=function(x,y,z,...)
  {
    panel.levelplot(x,y,z,...)
    ltext(x, y, labels=pvals, cex=labcex, font=1)
  }

  levelplot(df, panel=labels, xlab="", ylab="", pretty=TRUE, scales=list(x=list(rot=45, at=seq(1,ncol(df),1)), y=list(rot=180, at=seq(1,ncol(df),1))), aspect="fill", col.regions=cols, main=strTitle, cuts=100, at=seq(iLowerRange,iUpperRange,0.1))
}