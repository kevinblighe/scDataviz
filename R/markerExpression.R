markerExpression <- function(
  data,
  markers = NULL)
{
  if (!is.null(markers)) {
    mat <- data$expression[,which(colnames(data$expression) %in% markers)]
  } else {
    mat <- data$expression
  }

  for (i in 1:ncol(mat)) {      
    require(RColorBrewer)
    col <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
    col <- col(100)[as.numeric(cut(mat[,i], breaks = 1000))]
    sortidx <- order(mat[,i], decreasing = FALSE)
    plot(
      data$layout[sortidx,],
      main = paste(colnames(mat)[i], "expression"),
      xlab = "UMAP 1",
      ylab = "UMAP 2",
      xlim = c(min(data$layout[,1], na.rm = TRUE) - 1, max(data$layout[,1], na.rm = TRUE) + 1),
      ylim = c(min(data$layout[,2], na.rm = TRUE) - 1, max(data$layout[,2], na.rm = TRUE) + 1),
      col = col[sortidx],
      pch = ".")

    require(SDMTools)
    par(xpd = TRUE)
    points <- cbind(
      x = c(par("usr")[2] - 0.9, par("usr")[2] + 0.1, par("usr")[2] - 0.9, par("usr")[2] + 0.1),
      y = c(par("usr")[3], par("usr")[3], par("usr")[4], par("usr")[4]))
    legend.gradient(points,
      cols = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
      title = paste(""),
      limits = c(
        round(min(mat[,i], na.rm = TRUE), 2),
        round(max(mat[,i], na.rm = TRUE), 2)),
      cex = 1.2)
  }
}