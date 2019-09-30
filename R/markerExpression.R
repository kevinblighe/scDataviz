markerExpression <- function(
  data,
  marker = NULL)
{
  if (is.null(marker)) {
    for (i in 1:ncol(data$data)) {
      mat <- data$data
      require(RColorBrewer)
      col <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
      col <- col(100)[as.numeric(cut(mat[,i], breaks = 1000))]
      sortidx <- order(mat[,i], decreasing = FALSE)
      plot(
        data$layout[sortidx,],
        main = paste(colnames(data$data)[i], "expression"),
        xlab = "UMAP 1",
        ylab = "UMAP 2",
        xlim = c(-15, 10),
        ylim = c(-13, 14),
        col = col[sortidx],
        pch = ".")

      require(SDMTools)
      require(RColorBrewer)
      par(xpd = TRUE)
      points <- cbind(
        x = c(par("usr")[2] - 0.9, par("usr")[2] + 0.1, par("usr")[2] - 0.9, par("usr")[2] + 0.1),
        y = c(par("usr")[3], par("usr")[3], par("usr")[4], par("usr")[4]))
      legend.gradient(points,
        cols = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
        title = paste(""),
        limits = c(
          round(min(data$data[,i], na.rm = TRUE), 2),
          round(max(data$data[,i], na.rm = TRUE), 2)),
        cex = 1.2)
    }
  } else {
    i <- which(colnames(data$data) == marker)

    mat <- data$data
    require(RColorBrewer)
    col <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
    col <- col(100)[as.numeric(cut(mat[,i], breaks = 1000))]
    sortidx <- order(mat[,i], decreasing = FALSE)
    plot(
      data$layout[sortidx,],
      main = paste(colnames(data$data)[i], "expression"),
      xlab = "UMAP 1",
      ylab = "UMAP 2",
      xlim = c(-15, 10),
      ylim = c(-13, 14),
      col = col[sortidx],
      pch = ".")

    require(SDMTools)
    require(RColorBrewer)
    par(xpd = TRUE)
    points <- cbind(
      x = c(par("usr")[2] - 0.9, par("usr")[2] + 0.1, par("usr")[2] - 0.9, par("usr")[2] + 0.1),
      y = c(par("usr")[3], par("usr")[3], par("usr")[4], par("usr")[4]))
    legend.gradient(points,
      cols = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(1000),
      title = paste(""),
      limits = c(
        round(min(data$data[,i], na.rm = TRUE), 2),
        round(max(data$data[,i], na.rm = TRUE), 2)),
      cex = 1.2)
  }
}