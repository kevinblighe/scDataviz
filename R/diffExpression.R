diffExpression <- function(
  data,
  cluster1,
  cluster2,
  test = 'wilcox')
{
  c1 <- as.character(which(data$nnc$res.0.01 == cluster1))
  c2 <- as.character(which(data$nnc$res.0.01 == cluster2))
  mat <- t(data$data)
  colnames(mat) <- 1:ncol(mat)
  FindMarkers(
    object = mat,
    cells.1 = c1,
    cells.2 = c2,
    features = rownames(mat),
    #reduction = NULL,
    logfc.threshold = 0.0,
    test.use = test,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    only.pos = FALSE,
    max.cells.per.ident = Inf,
    random.seed = 1,
    latent.vars = NULL,
    min.cells.feature = 3,
    min.cells.group = 3,
    pseudocount.use = 1)
}