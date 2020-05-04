## ---- echo = FALSE, message = FALSE-------------------------------------------

  library(knitr)
  library(kableExtra)
  opts_chunk$set(tidy = FALSE, message = FALSE, warning = FALSE)


## ----getPackage, eval = FALSE-------------------------------------------------
#  
#      if (!requireNamespace('BiocManager', quietly = TRUE))
#          install.packages('BiocManager')
#  
#      BiocManager::install('scDataviz')
#  

## ----getPackageDevel, eval = FALSE--------------------------------------------
#  
#      devtools::install_github('kevinblighe/scDataviz')
#  

## ----Load, message = FALSE----------------------------------------------------

  library(scDataviz)


## ----readFCS, eval = FALSE----------------------------------------------------
#  
#    filelist <- list.files(
#      path = "scDataviz_data/FCS/",
#      pattern = "*.fcs|*.FCS",
#      full.names = TRUE)
#    filelist
#  
#    metadata <- data.frame(
#      sample = gsub('\\ [A-Za-z0-9]*\\.fcs$', '',
#        gsub('scDataviz_data\\/FCS\\/\\/', '', filelist)),
#      group = c(rep('Healthy', 7), rep('Disease', 11)),
#      treatment = gsub('\\.fcs$', '',
#        gsub('scDataviz_data\\/FCS\\/\\/[A-Z0-9]*\\ ', '', filelist)),
#      row.names = filelist,
#      stringsAsFactors = FALSE)
#    metadata
#  
#    sce <- processFCS(
#      files = filelist,
#      metadata = metadata,
#      transformation = TRUE,
#      transFun = function (x) asinh(x),
#      asinhFactor = 5,
#      downsample = 9000,
#      downsampleVar = 0.2,
#      newColnames = paste0('CD', 1:65))
#  

## ----load---------------------------------------------------------------------

  load(system.file('extdata/', 'complosome.rdata', package = 'scDataviz'))
  

## ----ex1, fig.height = 7, fig.width = 8, fig.cap = "Perform principal component analysis"----

  library(PCAtools)
  p <- pca(assay(sce, 'scaled'), metadata = metadata(sce))
  biplot(p,
    lab = NULL,
    pointSize = 0.5,
    colby = 'treatment',
    legendPosition = 'right',
    title = 'PCA applied to CyTOF data',
    caption = '9000 cells randomly selected after having filtered for low variance')


## ----addPCAdim----------------------------------------------------------------

  reducedDim(sce, 'PCA') <- p$rotated[,1:20]


## ----performUMAP--------------------------------------------------------------

  sce <- performUMAP(sce)


## ----elbowHorn----------------------------------------------------------------

  elbow <- findElbowPoint(p$variance)
  horn <- parallelPCA(assay(sce, 'scaled'))

  elbow
  horn$n


## ----performUMAP_PCA----------------------------------------------------------

  sce <- performUMAP(sce, reducedDim = 'PCA', dims = c(1:horn$n))


## ----ex2, fig.height = 8, fig.width = 16, fig.cap = "Create a contour plot of the UMAP layout"----

  ggout1 <- contourPlot(sce,
    reducedDim = 'UMAP',
    subtitle = 'UMAP performed on expression values',
    legendLabSize = 18,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)

  ggout2 <- contourPlot(sce,
    reducedDim = 'UMAP_PCA',
    subtitle = 'UMAP performed on PC eigenvectors',
    legendLabSize = 18,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)

  plot_grid(ggout1, ggout2,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 24)


## ----ex3, fig.height = 12, fig.width = 20, fig.cap = "Show marker expression across the layout"----

  markers <- sample(rownames(sce), 6)
  markers

  ggout1 <- markerExpression(sce,
    markers = markers,
    subtitle = 'UMAP performed on expression values',
    nrow = 1, ncol = 6,
    legendKeyHeight = 1.0,
    legendLabSize = 18,
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)

  ggout2 <-  markerExpression(sce,
    markers = markers,
    reducedDim = 'UMAP_PCA',
    subtitle = 'UMAP performed on PC eigenvectors',
    nrow = 1, ncol = 6,
    legendKeyHeight = 1.0,
    legendLabSize = 18,
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)

  plot_grid(ggout1, ggout2,
    labels = c('A','B'),
    nrow = 2, align = "l", label_size = 24)


## ----metadataPlot-------------------------------------------------------------

  head(metadata(sce))

  levels(metadata(sce)$group)

  levels(metadata(sce)$treatment)


## ----ex4, fig.height = 12, fig.width = 14, fig.cap = "Shade cells by metadata", message = FALSE----

  ggout1 <- metadataPlot(sce,
    colby = 'group',
    colkey = c(Healthy = 'royalblue', Disease = 'red2'),
    title = 'Disease status',
    subtitle = 'UMAP performed on expression values',
    legendLabSize = 16,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  ggout2 <- metadataPlot(sce,
    reducedDim = 'UMAP_PCA',
    colby = 'group',
    colkey = c(Healthy = 'royalblue', Disease = 'red2'),
    title = 'Disease status',
    subtitle = 'UMAP performed on PC eigenvectors',
    legendLabSize = 16,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  ggout3 <- metadataPlot(sce,
    colby = 'treatment',
    title = 'Treatment type',
    subtitle = 'UMAP performed on expression values',
    legendLabSize = 16,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  ggout4 <- metadataPlot(sce,
    reducedDim = 'UMAP_PCA',
    colby = 'treatment',
    title = 'Treatment type',
    subtitle = 'UMAP performed on PC eigenvectors',
    legendLabSize = 16,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  plot_grid(ggout1, ggout3, ggout2, ggout4,
    labels = c('A','B','C','D'),
    nrow = 2, ncol = 2, align = "l", label_size = 24)


## ----ex5, message = FALSE, fig.height = 8, fig.width = 14, fig.cap = "Find ideal clusters in the UMAP layout via k-nearest neighbours"----

  sce <- clusKNN(sce,
    k.param = 20,
    prune.SNN = 1/15,
    resolution = 0.01,
    algorithm = 2)

  sce <- clusKNN(sce,
    reducedDim = 'UMAP_PCA',
    clusterAssignName = 'Cluster_PCA',
    k.param = 20,
    prune.SNN = 1/15,
    resolution = 0.01,
    algorithm = 2)

  ggout1 <- plotClusters(sce,
    clusterColname = 'Cluster',
    labSize = 7.0,
    subtitle = 'UMAP performed on expression values',
    caption = paste0('Note: clusters / communities identified via',
      '\nLouvain algorithm with multilevel refinement'),
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  ggout2 <- plotClusters(sce,
    clusterColname = 'Cluster_PCA',
    reducedDim = 'UMAP_PCA',
    labSize = 7.0,
    subtitle = 'UMAP performed on PC eigenvectors',
    caption = paste0('Note: clusters / communities identified via',
      '\nLouvain algorithm with multilevel refinement'),
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

  plot_grid(ggout1, ggout2,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 24)


## ----ex6a, fig.height = 9, fig.width = 18, fig.cap = "Plot marker expression per identified cluster1"----

  markerExpressionPerCluster(sce,
    caption = 'Cluster assignments based on UMAP performed on expression values',
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)


## ----ex6b, fig.height = 9, fig.width = 16, fig.cap = "Plot marker expression per identified cluster2"----

  clusters <- unique(metadata(sce)[['Cluster_PCA']])
  clusters

  markers <- sample(rownames(sce), 8)
  markers

  markerExpressionPerCluster(sce,
    clusters = clusters,
    clusterAssign = metadata(sce)[['Cluster_PCA']],
    markers = markers,
    nrow = 2, ncol = 4,
    caption = 'Cluster assignments based on UMAP performed on PC eigenvectors',
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)


## ----ex7, fig.height = 8, fig.width = 16, fig.cap = "Determine enriched markers in each cluster and plot the expression signature"----

  plotSignatures(sce,
    labCex = 1.5,
    legendCex = 1.5)


## ----importRandomData1--------------------------------------------------------

  mat <- jitter(matrix(
    MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
    ncol = 20))
  colnames(mat) <- paste0('CD', 1:ncol(mat))
  rownames(mat) <- paste0('cell', 1:nrow(mat))

  metadata <- data.frame(
    group = rep('A', nrow(mat)),
    row.names = rownames(mat),
    stringsAsFactors = FALSE)
  head(metadata)

  sce <- importData(mat,
    assayname = 'normcounts',
    metadata = metadata)
  sce


## ----importRandomData2--------------------------------------------------------

  sce <- importData(mat,
    assayname = 'normcounts',
    metadata = NULL)
  sce


## -----------------------------------------------------------------------------

sessionInfo()


