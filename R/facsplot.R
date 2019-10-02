facsplot <- function(
  data,
  marker1,
  marker2,
  log2 = TRUE,
  bins = 200,
  main = 'FACS-like plot',
  x1 = 0,
  x2 = 0,
  y1 = 0,
  y2 = 0,
  cex = 1,
  colramp = colorRampPalette(rev(brewer.pal(9,"Spectral"))))
{
  data <- as.data.frame(data$expression)

  x <- data[[marker1]]
  y <- data[[marker2]]

  x[x==0] <- NA
  if (log2) {x <- log2(x)}

  y[y==0] <- NA
  if (log2) {y <- log2(y)}

  xlim <- c(min(x), max(x))
  ylim <- c(min(y), max(y))

  hbplot <- hexbinplot(
    y ~ x,
    data = data,
    aspect = '1',
    xbins = bins,
    xlab = marker1,
    ylab = marker2,
    xlim = xlim,
    ylim = ylim,
    cex.labels = cex,
    cex.title = cex,
    colramp = colramp,
    panel=function(x, y, ...) {
      panel.hexbinplot(x, y, ...)
      panel.abline(v=c(x1,x2), h=c(y1,y2), col="black", lwd=4, lty=5)})

  return(hbplot)
}