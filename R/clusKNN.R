clusKNN <- function(
  sce,
  reducedDim = 'UMAP',
  dimColnames = c('UMAP1','UMAP2'),
  clusterAssignName = 'Cluster',
  distance.matrix = FALSE,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "rann",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  modularity.fxn = 1,
  initial.membership = NULL,
  weights = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  method = "matrix",
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL)
{
  layout <- reducedDim(sce, reducedDim)[,dimColnames]

  layout$nn <- FindNeighbors(
    layout,
    distance.matrix = distance.matrix,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    force.recalc = force.recalc)

  layout$nnc <- FindClusters(
    layout$nn$snn,
    modularity.fxn = modularity.fxn,
    initial.membership = initial.membership,
    weights = node.sizes,
    node.sizes = node.sizes,
    resolution = resolution,
    method = method,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    group.singletons = group.singletons,
    temp.file.location = temp.file.location,
    edge.file.name = edge.file.name,
    verbose = verbose)

  if (length(which(colnames(metadata(sce)) == clusterAssignName)) > 0 ) {
    metadata(sce) <- metadata(sce)[,-which(colnames(metadata(sce)) == clusterAssignName)]
  }

  metadata(sce) <- data.frame(
    metadata(sce),
    as.numeric(as.character(layout$nnc[,1])))

  colnames(metadata(sce))[ncol(metadata(sce))] <- clusterAssignName

  return(sce)
}

