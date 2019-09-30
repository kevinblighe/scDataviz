clusKNN <- function(
  data)
{
  require(Seurat)

  data$nn <- FindNeighbors(
    data$layout,
    distance.matrix = FALSE,
    k.param = 10,
    compute.SNN = TRUE,
    prune.SNN = 1/15,
    nn.eps = 0,
    verbose = TRUE,
    force.recalc = FALSE)
  # @param modularity.fxn Modularity function (1 = standard; 2 = alternative).
  # @param initial.membership,weights,node.sizes Parameters to pass to the Python leidenalg function.
  # @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
  # @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
  # @param n.start Number of random starts.
  # @param n.iter Maximal number of iterations per random start.
  # @param random.seed Seed of the random number generator.
  # @param temp.file.location Directory where intermediate files will be written. Specify the ABSOLUTE path.
  # @param edge.file.name Edge file to use as input for modularity optimizer jar.
  # @param verbose Print output
  data$nnc <- FindClusters(
    data$nn$snn,
    modularity.fxn = 1,
    initial.membership = NULL,
    weights = NULL,
    node.sizes = NULL,
    resolution = 0.01,
    algorithm = 2, #1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    temp.file.location = NULL,
    edge.file.name = NULL,
    verbose = TRUE)

  data$colour <- colorRampPalette(brewer.pal(9,"Spectral"))(length(unique(data$nnc[,1])))[data$nnc[,1]]
  data$lab <- as.character(data$nnc[,1])
  data$lab[duplicated(data$lab)] <- ""

  return(data)
}
