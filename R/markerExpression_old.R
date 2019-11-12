markerExpression <- function(
  sce,
  markers = NULL)
{
  data <- as.data.frame(t(assay(sce, 'scaled')))

  u <- as.data.frame(reducedDim(sce, "UMAP"))

  if (!is.null(markers)) {
    data <- data[,which(colnames(data) %in% markers)]
  } else {
    data <- data
  }

  for (i in 1:ncol(data)) {
    col <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
    col <- col(100)[as.numeric(cut(data[,i], breaks = 1000))]
    sortidx <- order(data[,i], decreasing = FALSE)
    plot(
      u[sortidx,],
      main = paste(colnames(data)[i], "expression"),
      xlab = "UMAP 1",
      ylab = "UMAP 2",
      xlim = c(min(u[,1], na.rm = TRUE) - 1, max(u[,1], na.rm = TRUE) + 1),
      ylim = c(min(u[,2], na.rm = TRUE) - 1, max(u[,2], na.rm = TRUE) + 1),
      col = col[sortidx],
      pch = ".")

    par(xpd = TRUE)
    points <- cbind(
      x = c(par("usr")[2] - 0.9, par("usr")[2] + 0.1, par("usr")[2] - 0.9, par("usr")[2] + 0.1),
      y = c(par("usr")[3], par("usr")[3], par("usr")[4], par("usr")[4]))
    legend.gradient(points,
      cols = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
      title = paste(""),
      limits = c(
        round(min(data[,i], na.rm = TRUE), 2),
        round(max(data[,i], na.rm = TRUE), 2)),
      cex = 1.2)
  }
}
