clusKNN <- function(
  indata,
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
  edge.file.name = NULL,
  overwrite = FALSE)
{

  if (is(indata, 'SingleCellExperiment')) {
    if (verbose) message('--input data class is SingleCellExperiment')
    if (verbose) message('--\'', reducedDim, '\' reduced dimensional component will ',
      'be used for clustering, with dims / columns: ', paste(dimColnames,
        collapse = ', '))
    clusdata <- reducedDim(indata, reducedDim)[,dimColnames]

    if (length(which(colnames(metadata(indata)) == clusterAssignName)) && overwrite == FALSE) {
      stop('Column \'', clusterAssignName, '\' already found in metadata.. ',
        'Re-run the command with \'overwrite == TRUE\' in order to confirm ',
        'overwrite, or change the value of \'clusterAssignName\'.')
    }
  } else {
    if (verbose) message('--input data class is ', class(indata))
    clusdata <- indata
  }

  nn <- FindNeighbors(
    clusdata,
    distance.matrix = distance.matrix,
    k.param = k.param,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    force.recalc = force.recalc)

  nnc <- FindClusters(
    nn$snn,
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

  if (is(indata, 'SingleCellExperiment')) {
    if (length(which(colnames(metadata(indata)) == clusterAssignName))) {
      metadata(indata) <- metadata(indata)[,-which(colnames(metadata(indata)) == clusterAssignName)]
    }

    if (length(metadata(indata))) {
      metadata(indata) <- data.frame(
        metadata(indata),
        as.numeric(as.character(nnc[,1])),
        row.names = rownames(metadata(indata)))
    } else {
      metadata(indata) <- data.frame(
        as.numeric(as.character(nnc[,1])),
        row.names = rownames(metadata(indata)))
    }

    colnames(metadata(indata))[ncol(metadata(indata))] <- clusterAssignName

    if (verbose) message('\ncluster information added to your input SingleCellExperiment ',
      'object\'s metadata under colname \'', clusterAssignName, '\'')

    return(indata)

  } else {

    return(as.numeric(as.character(nnc[,1])))

  }
}
