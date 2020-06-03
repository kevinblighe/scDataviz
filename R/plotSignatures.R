#' Find enriched markers per identified cluster and visualise these as a custom corrplot.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'    data-frame or matrix, this should relate to expression data (cells as
#'    columns; genes as rows). If a \code{SingleCellExperiment} object, data will be
#'    extracted from an assay component named by \code{assay}.
#' @param assay Name of the assay slot in \code{indata} from which data will be
#'    taken, assuming \code{indata} is a \code{SingleCellExperiment} object.
#' @param clusterAssign A vector of cell-to-cluster assignments. This can be
#'    from any source but must align with your cells / variables. There is no
#'    check to ensure this when \code{indata} is not a \code{SingleCellExperiment} object.
#' @param funcSummarise A mathematical function used to summarise expression
#'    per marker, per cluster.
#' @param col colorRampPalette to be used for shading low-to-high expression.
#' @param labCex cex (size) of the main plot labels.
#' @param legendCex cex (size) of the legend labels.
#' @param labDegree Rotation angle of the main plot labels.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Find enriched markers per identified cluster and visualise these as a custom corrplot. \code{plotSignatures} first collapses your input data's expression profiles from the level of cells to the level of clusters based on a mathematical function specified by \code{funcSummarise}. It then centers and scales the data range to be between -1 and +1 for visualisation purposes.
#'
#' @return A \code{corrplot} object.
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
#' u <- umap::umap(mat)$layout
#' colnames(u) <- c('UMAP1','UMAP2')
#' rownames(u) <- rownames(mat)
#' clus <- clusKNN(u)
#'
#' plotSignatures(t(mat), clusterAssign = clus)
#'
#' @import SingleCellExperiment
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom stats aggregate median sd
#' @importFrom corrplot corrplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend
#' @importFrom scales rescale
#' @importFrom methods is
#'
#' @export
plotSignatures <- function(
  indata,
  assay = 'scaled',
  clusterAssign = metadata(indata)[['Cluster']],
  funcSummarise = function(x) mean(x, na.rm = TRUE),
  col = colorRampPalette(brewer.pal(9, 'RdPu'))(100),
  labCex = 1.0,
  legendCex = 1.0,
  labDegree = 90,
  verbose = TRUE)
{

  if (is(indata, 'SingleCellExperiment')) {
    if (verbose) message('--input data class is SingleCellExperiment')
    data <- as.data.frame(t(assay(indata, assay)))
  } else {
    if (verbose) message('--input data class is ', class(indata))
    data <- as.data.frame(t(indata))
  }

  data <- aggregate(data, list(clusterAssign), funcSummarise)
  rownames(data) <- data[,1]
  data <- data.matrix(data[,-1])
  data <- t(scale(t(data), center = TRUE, scale = FALSE))
  data <- apply(data, 1, function(x) rescale(x, c(-1,1)))

  corrplot(
    data.matrix(t(data)),
    method = 'circle',
    order = 'original',
    addgrid.col = 'grey60',
    tl.col = 'black',
    col = col,
    cl.cex = labCex,
    cl.pos = 'b',
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
