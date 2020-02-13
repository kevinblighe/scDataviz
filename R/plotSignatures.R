plotSignatures <- function(
  indata,
  assay = 'scaled',
  clusterAssign = metadata(indata)[['Cluster']],
  funcSummarise = function(x) median(x, na.rm = TRUE),
  col = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100),
  labCex = 1.0,
  legendCex = 1.0,
  labDegree = 90)
{
  if (class(indata) == 'SingleCellExperiment') {
    message('--input data class is SingleCellExperiment')
    data <- as.data.frame(t(as.matrix(assay(indata, assay))))
  } else {
    message('--input data class is ', class(indata))
    data <- as.data.frame(t(as.matrix(indata)))
  }

  data <- aggregate(data, list(clusterAssign), funcSummarise)
  data <- data[,-1]
  data <- apply(data, 2, scale, scale = FALSE)
  #data <- t(data) / max(abs(range(data)))
  data <- rescale(data, c(-1,1))

  corrplot(
    data.matrix(data),
    method = "circle",
    order = "original",
    addgrid.col = "grey60",
    tl.col = "black",
    col = col,
    cl.cex = labCex,
    cl.pos = "b",
    cl.ratio = 0.4,
    cl.lim = c(-1,1),
    cl.length = 3,
    tl.cex = labCex,
    tl.srt = labDegree,
    mar = c(1,3,1,2))
  legend(
    'top',
    bty = 'n',
    cex = legendCex,
    title = '',
    c('High', 'Average', 'Low'),
    fill = c(col[length(col)-10], col[length(col)/2], col[11]),
    horiz = TRUE)
}
