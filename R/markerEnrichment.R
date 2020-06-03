#' Find enriched markers per identified cluster and calculate cluster abundances across these for samples and metadata variables.
#'
#' @param indata A data-frame or matrix, or \code{SingleCellExperiment} object. If a
#'   data-frame or matrix, this should relate to expression data (cells as
#'   columns; genes as rows). If a \code{SingleCellExperiment} object, data will be
#'   extracted from an assay component named by \code{assay}.
#' @param meta If 'indata' is a non-\code{SingleCellExperiment} object, \code{meta} must be
#'   activated and relate to a data-frame of metadata that aligns with the columns
#'   of \code{indata}, and that also contains a column name specified by \code{studyvarID}.
#' @param assay Name of the assay slot in \code{indata} from which data will be
#'   taken, assuming \code{indata} is a \code{SingleCellExperiment} object.
#' @param sampleAbundances Logical, indicating whether or not to calculate
#'   cluster abundances across study samples.
#' @param sampleID If \code{sampleAbundances == TRUE}, a column name from the
#'   provided metadata representing over which sample cluster abundances
#'   will be calculated.
#' @param studyvarID A column name from the provided metadata representing a
#'   condition or trait over which cluster abundances will be calculated.
#' @param clusterAssign A vector of cell-to-cluster assignments. This can be
#'   from any source but must align with your cells / variables. There is no
#'   check to ensure this when 'indata' is not a \code{SingleCellExperiment}
#'  object.
#' @param funcSummarise A mathematical function used to summarise expression
#'   per marker per cluster.
#' @param method Type of summarisation to apply to the data for final marker
#'   selection. Possible values include \code{Z} or \code{quantile}. If
#'   \code{Z}, \code{limits} relate to lower and upper Z-score cut-offs for
#'   low|high markers. The defaults of -1.96 and +1.96 are equivalents of
#'   p<0.05 on a two-tailed distribution. If \code{quantile}, \code{prob}
#'   will be used to define the \code{n}th lower and 1 - \code{n}th upper
#'   quantiles, which will be used for selecting low|high markers.
#' @param prob See details for \code{method}.
#' @param limits See details for \code{method}.
#' @param verbose Boolean (TRUE / FALSE) to print messages to console or not.
#'
#' @details
#' Find enriched markers per identified cluster and calculate cluster abundances across these for samples and metadata variables. \code{markerEnrichment} first collapses your input data's expression profiles from the level of cells to the level of clusters based on a mathematical function specified by \code{funcSummarise}. It then either selects, per cluster, low|high markers via quantiles, or transforms this collapsed data to global Z-scores and selects low|high markers based on Z-score cut-offs.
#'
#' @return A \code{data.frame} object.
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
#' metadata <- data.frame(
#'   group = c(rep('PB1', 25), rep('PB2', 25)),
#'   row.names = rownames(u))
#'
#' markerEnrichment(t(mat), meta = metadata,
#'   sampleAbundances = FALSE,
#'   studyvarID = 'group', clusterAssign = clus)
#'
#' @import SingleCellExperiment
#'
#' @importFrom MASS rnegbin
#' @importFrom umap umap
#' @importFrom stats aggregate median sd quantile
#' 
#' @export
markerEnrichment <- function(
  indata,
  meta = NULL,
  assay = 'scaled',
  sampleAbundances = TRUE,
  sampleID = 'sample',
  studyvarID = NULL,
  clusterAssign = metadata(indata)[['Cluster']],
  funcSummarise = function(x) mean(x, na.rm = TRUE),
  method = 'Z',
  prob = 0.1,
  limits = c(-1.96, 1.96),
  verbose = TRUE)
{
  iCellsPerCluster <- c()
  iTotalCells <- c()
  iPercentage <- c()
  NegativeMarkers <- c()
  PositiveMarkers <- c()

  if (is(indata, 'SingleCellExperiment')) {

    if (verbose) message('--input data class is SingleCellExperiment')
    metadata <- metadata(indata)
    data <- as.data.frame(t(as.matrix(assay(indata, assay))))

  } else {

    if (verbose) message('--input data class is ', class(indata))

    if (is.null(meta)) {
      stop('When the input data is a non-SingleCellExperiment object, ',
        '\'indata\' must relate to an expression matrix (cells as columns; ',
        'genes as rows), while \'meta\' must be non-NULL and relate to ',
        'metadata assocaited with this data.')
    } else if (!all(rownames(meta) == colnames(indata))) {
      stop('\'rownames(meta)\' must be equal to \'colnames(indata)\'')
    }

    metadata <- meta
    data <- as.data.frame(t(as.matrix(indata)))
  }

  if (sampleAbundances) if (verbose)
    message('--sample cluster abundances will be determined based on \'',
      sampleID, '\' metadata column')

  if (!is.null(studyvarID)) if (verbose)
    message('--cluster abundances will be determined for \'',
      studyvarID, '\' variable from metadata')

  data <- aggregate(data, list(clusterAssign), funcSummarise)
  data <- data.matrix(data[,-1])

  if (method == 'quantile') {
    if (verbose) message('--marker selection method is \'quantile\'')
    quarts <- quantile(data, probs = seq(0, 1, prob))
    limits[1] <- quarts[2]
    limits[2] <- quarts[length(quarts) -1]
    data <- t(data)
  } else if (method == 'Z') {
    if (verbose) message('--marker selection method is \'Z\'')
    data <- t((data - mean(data, na.rm = TRUE)) / sd(data, na.rm = TRUE))
  }

  nclus <- length(unique(clusterAssign))

  res <- data.frame(row.names = seq(0, nclus-1))
  metares <- data.frame(row.names = seq(0, nclus-1))

  for (j in seq(0, nclus-1)) {
    iCellsPerCluster <- length(clusterAssign[clusterAssign == j])
    iTotalCells <- length(clusterAssign)
    iPercentage <- (iCellsPerCluster/iTotalCells) * 100

    NegativeMarkers <- names(which(data[,j+1] <= limits[1]))
    PositiveMarkers <- names(which(data[,j+1] >= limits[2]))

    if (sampleAbundances) {
      metaclusAbundance <- table(metadata[which(clusterAssign == j), sampleID])
      metaclusAbundance <- (metaclusAbundance / iCellsPerCluster) * 100
      metaclusAbundance <- metaclusAbundance[match(sort(unique(metadata[,sampleID])), names(metaclusAbundance))]
    }

    if (!is.null(studyvarID)) {
        studyvar <- table(metadata[which(clusterAssign == j),studyvarID])
    }

    res <- rbind(res,
      c(j,
        iCellsPerCluster,
        iTotalCells,
        iPercentage,
        paste(NegativeMarkers, '-', sep = '', collapse = ''),
        paste(PositiveMarkers, '+', sep = '', collapse = '')))
    colnames(res) <- c('Cluster', 'nCells','TotalCells','PercentCells','NegMarkers','PosMarkers')
    res$Cluster <- as.numeric(as.character(res$Cluster))
    res$nCells <- as.numeric(as.character(res$nCells))
    res$TotalCells <- as.numeric(as.character(res$TotalCells))
    res$PercentCells <- as.numeric(as.character(res$PercentCells))
    res$PosMarkers <- as.character(res$PosMarkers)
    res$NegMarkers <- as.character(res$NegMarkers)

    if (sampleAbundances & !is.null(studyvarID)) {
      metares <- rbind(metares,
        c(metaclusAbundance, studyvar))
      colnames(metares) <- c(paste0('PerCent_', names(metaclusAbundance)), paste0('nCell_', names(studyvar)))
    } else if (sampleAbundances & is.null(studyvarID))  {
      metares <- rbind(metares,
        c(metaclusAbundance))
      colnames(metares) <- c(paste0('PerCent_', names(metaclusAbundance)))
    } else if (!sampleAbundances & !is.null(studyvarID))  {
      metares <- rbind(metares,
        c(studyvar))
      colnames(metares) <- c(paste0('nCell_', names(studyvar)))
    }
  }

  res$PosMarkers[res$PosMarkers == '+'] <- NA
  res$NegMarkers[res$NegMarkers == '-'] <- NA

  return(cbind(res, metares))
}
