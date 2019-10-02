plotSignatures <- function(
  data,
  col = colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100),
  cexlab = 1.0,
  cexlegend = 1.0,
  labDegree = 90)
{
  df <- data$expression
  df <- aggregate(data.matrix(data$expression), data$nnc, mean)
  df <- df[,-1]
  df <- df[,order(colnames(df))]

  # center the medoids
  df <- apply(df, 2, scale, scale = FALSE)

  # set max to +1
  df <- t(df) / max(abs(range(df)))

  # stretch out to -1 / +1 scale
  df <- scales::rescale(df, c(-1,1))

  corrplot(
    data.matrix(t(df)),
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