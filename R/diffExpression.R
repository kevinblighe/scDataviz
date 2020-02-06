diffExpression <- function(
  sce,
  cells1,
  cells2,
  features = rownames(sce),
  logfc.threshold = 0.0,
  test.use = 'wilcox',
  min.pct = 0.0,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1)
{
  data <- as.data.frame(assay(sce, 'scaled'))

  FindMarkers(
    object = data,
    cells.1 = cells1,
    cells.2 = cells2,
    features = features,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    latent.vars = latent.vars,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    pseudocount.use = pseudocount.use)
}
