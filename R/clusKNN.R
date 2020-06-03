#' A wrapper function for \code{Seurat}'s \code{FindNeighbors} and \code{FindClusters}.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object.
#'   If a \code{SingleCellExperiment} object, the cell-to-cluster assignments
#'   will be added as a new column, specified by \code{clusterAssignName}, to
#'   the input object's metadata; if a data-frame or matrix, only the cluster
#'   assignment vector is returned.
#' @param reducedDim A reduced dimensional component stored within \code{indata}.
#'   e.g., PCA or UMAP.
#' @param dimColnames The column names of the dimensions to use.
#' @param clusterAssignName The new column name in the metadata that will
#'   contain the determined cell-to-cluster assignments.
#' @param distance.matrix Refer to \code{?Seurat::FindNeighbors}.
#' @param k.param Refer to \code{?Seurat::FindNeighbors}.
#' @param compute.SNN Refer to \code{?Seurat::FindNeighbors}.
#' @param prune.SNN Refer to \code{?Seurat::FindNeighbors}.
#' @param nn.method Refer to \code{?Seurat::FindNeighbors}.
#' @param annoy.metric Refer to \code{?Seurat::FindNeighbors}.
#' @param nn.eps Refer to \code{?Seurat::FindNeighbors}.
#' @param verbose Refer to \code{?Seurat::FindNeighbors}.
#' @param force.recalc Refer to \code{?Seurat::FindNeighbors}.
#' @param modularity.fxn Refer to \code{?Seurat::FindClusters}.
#' @param initial.membership Refer to \code{?Seurat::FindClusters}.
#' @param weights Refer to \code{?Seurat::FindClusters}.
#' @param node.sizes Refer to \code{?Seurat::FindClusters}.
#' @param resolution Refer to \code{?Seurat::FindClusters}.
#' @param method Refer to \code{?Seurat::FindClusters}.
#' @param algorithm Refer to \code{?Seurat::FindClusters}.
#' @param n.start Refer to \code{?Seurat::FindClusters}.
#' @param n.iter Refer to \code{?Seurat::FindClusters}.
#' @param random.seed Refer to \code{?Seurat::FindClusters}.
#' @param group.singletons Refer to \code{?Seurat::FindClusters}.
#' @param temp.file.location Refer to \code{?Seurat::FindClusters}.
#' @param edge.file.name Refer to \code{?Seurat::FindClusters}.
#' @param overwrite When the input object is a SingleCellExperiment, enabling
#'   this will result in the overwriting, with the new cluster assignments, of
#'   any column in your metadata that has the same name as
#'   \code{clusterAssignName}.
#'
#' @details
#' A wrapper function for Seurat's FindNeighbors and FindClusters.
#'
#' @return A \code{SingleCellExperiment} or \code{numeric} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create random data that follows a negative binomial
#' mat <- jitter(matrix(
#'   MASS::rnegbin(rexp(1000, rate=.1), theta = 4.5),
#'   ncol = 20))
#' colnames(mat) <- paste0('CD', 1:ncol(mat))
#' rownames(mat) <- paste0('cell', 1:nrow(mat))
#'
#' clusKNN(mat)
#'
#' @import SingleCellExperiment
#'
#' @importFrom MASS rnegbin
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' 
#' @export
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
