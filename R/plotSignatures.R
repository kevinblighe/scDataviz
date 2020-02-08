plotSignatures <- function(
  sce,
  assay = 'scaled',
  clusterVector = metadata(sce)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  col = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100),
  cexlab = 1.0,
  cexlegend = 1.0,
  labDegree = 90)
{
  data <- as.data.frame(t(assay(sce, assay)))
  data <- aggregate(data, list(clusterVector), funcSummarise)
  data <- data[,-1]
  data <- apply(data, 2, scale, scale = FALSE)
  data <- t(data) / max(abs(range(data)))
  data <- rescale(data, c(-1,1))

  corrplot(
    data.matrix(t(data)),
    method = "circle",
    order = "original",
    addgrid.col = "grey60",
    tl.col = "black",
    col = col,
    cl.cex = cexlab,
    cl.pos = "b",
    cl.ratio = 0.4,
    cl.lim = c(-1,1),
    cl.length = 3,
    tl.cex = cexlab,
    tl.srt = labDegree,
    mar = c(1,3,1,2))
  legend(
    'top',
    bty = 'n',
    cex = cexlegend,
    title = '',
    c('High', 'Average', 'Low'),
    fill = c(col[length(col)-10], col[length(col)/2], col[11]),
    horiz = TRUE)
}
