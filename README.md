scDataviz: single cell dataviz and downstream analyses
================
Kevin Blighe
2020-06-03

# Introduction

In the single cell World, which includes flow cytometry, mass cytometry,
single-cell RNA-seq (scRNA-seq), and others, there is a need to improve
data visualisation and to bring analysis capabilities to researchers
even from non-technical backgrounds. *scDataviz* (Blighe 2020) attempts
to fit into this space, while also catering for advanced users.
Additonally, due to the way that *scDataviz* is designed, which is based
on *SingleCellExperiment* (Lun and Risso 2020), it has a ‘plug and play’
feel, and immediately lends itself as flexibile and compatibile with
studies that go beyond *scDataviz*. Finally, the graphics in *scDataviz*
are generated via the *ggplot* (Wickham 2016) engine, which means that
users can ‘add on’ features to these with ease.

# Installation

## 1\. Download the package from Bioconductor

``` r
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

  BiocManager::install('scDataviz')
```

Note: to install development version:

``` r
  devtools::install_github('kevinblighe/scDataviz')
```

## 2\. Load the package into R session

``` r
  library(scDataviz)
```

# Tutorial 1: CyTOF FCS data

Here, we will start with sample data stored as FCS files. Specifically,
we will utilise some of the data used in [Deep phenotyping detects a
pathological CD4+ T-cell complosome signature in systemic
sclerosis](https://www.nature.com/articles/s41423-019-0360-8).

In order to download this data, we can use `git clone` from your command
prompt:

``` bash

  git clone https://github.com/kevinblighe/scDataviz_data ;
```

**NB - this command (above) needs to be run outside R at your shell’s
command prompt (e.g., BASH)**

Now, read in the data and normalise it. The `processFCS` function, by
default, removes variables based on low variance and also downsamples
\[randomly\] your data to 100000 variables. The user can change these
via the `downsample` and `downsampleVar` parameters.

``` r
  filelist <- list.files(
    path = "scDataviz_data/FCS/",
    pattern = "*.fcs|*.FCS",
    full.names = TRUE)
  filelist
```

    ##  [1] "scDataviz_data/FCS//HD00 CD46.fcs"   
    ##  [2] "scDataviz_data/FCS//HD00 Unstim.fcs" 
    ##  [3] "scDataviz_data/FCS//HD01 CD46.fcs"   
    ##  [4] "scDataviz_data/FCS//HD01 Unstim.fcs" 
    ##  [5] "scDataviz_data/FCS//HD262 CD3.fcs"   
    ##  [6] "scDataviz_data/FCS//HD262 CD46.fcs"  
    ##  [7] "scDataviz_data/FCS//HD262 Unstim.fcs"
    ##  [8] "scDataviz_data/FCS//P00 CD46.fcs"    
    ##  [9] "scDataviz_data/FCS//P00 Unstim.fcs"  
    ## [10] "scDataviz_data/FCS//P02 CD3.fcs"     
    ## [11] "scDataviz_data/FCS//P02 CD46.fcs"    
    ## [12] "scDataviz_data/FCS//P03 CD3.fcs"     
    ## [13] "scDataviz_data/FCS//P03 CD46.fcs"    
    ## [14] "scDataviz_data/FCS//P04 CD3.fcs"     
    ## [15] "scDataviz_data/FCS//P04 CD46.fcs"    
    ## [16] "scDataviz_data/FCS//P08 CD3.fcs"     
    ## [17] "scDataviz_data/FCS//P08 CD46.fcs"    
    ## [18] "scDataviz_data/FCS//P08 Unstim.fcs"

``` r
  metadata <- data.frame(
    sample = gsub('\\ [A-Za-z0-9]*\\.fcs$', '',
      gsub('scDataviz_data\\/FCS\\/\\/', '', filelist)),
    group = c(rep('Healthy', 7), rep('Disease', 11)),
    treatment = gsub('\\.fcs$', '',
      gsub('scDataviz_data\\/FCS\\/\\/[A-Z0-9]*\\ ', '', filelist)),
    row.names = filelist,
    stringsAsFactors = FALSE)
  metadata
```

    ##                                      sample   group treatment
    ## scDataviz_data/FCS//HD00 CD46.fcs      HD00 Healthy      CD46
    ## scDataviz_data/FCS//HD00 Unstim.fcs    HD00 Healthy    Unstim
    ## scDataviz_data/FCS//HD01 CD46.fcs      HD01 Healthy      CD46
    ## scDataviz_data/FCS//HD01 Unstim.fcs    HD01 Healthy    Unstim
    ## scDataviz_data/FCS//HD262 CD3.fcs     HD262 Healthy       CD3
    ## scDataviz_data/FCS//HD262 CD46.fcs    HD262 Healthy      CD46
    ## scDataviz_data/FCS//HD262 Unstim.fcs  HD262 Healthy    Unstim
    ## scDataviz_data/FCS//P00 CD46.fcs        P00 Disease      CD46
    ## scDataviz_data/FCS//P00 Unstim.fcs      P00 Disease    Unstim
    ## scDataviz_data/FCS//P02 CD3.fcs         P02 Disease       CD3
    ## scDataviz_data/FCS//P02 CD46.fcs        P02 Disease      CD46
    ## scDataviz_data/FCS//P03 CD3.fcs         P03 Disease       CD3
    ## scDataviz_data/FCS//P03 CD46.fcs        P03 Disease      CD46
    ## scDataviz_data/FCS//P04 CD3.fcs         P04 Disease       CD3
    ## scDataviz_data/FCS//P04 CD46.fcs        P04 Disease      CD46
    ## scDataviz_data/FCS//P08 CD3.fcs         P08 Disease       CD3
    ## scDataviz_data/FCS//P08 CD46.fcs        P08 Disease      CD46
    ## scDataviz_data/FCS//P08 Unstim.fcs      P08 Disease    Unstim

``` r
  sce <- processFCS(
    files = filelist,
    metadata = metadata,
    transformation = TRUE,
    transFun = function (x) asinh(x),
    asinhFactor = 5,
    downsample = 25000,
    downsampleVar = 0.2,
    newColnames = paste0('CD', 1:65))
```

One can also create a new *SingleCellExperiment* object manually using
any type of data, including any data from scRNA-seq produced elsewhere.
Import functions for data deriving from other sources is covered in
Tutorials 2 and 3 in this vignette. All functions in *scDataviz*
additionally accept data-frames or matrices on their own,
de-necessitating the reliance on the *SingleCellExperiment* class.

## Perform principal component analysis (PCA)

We can use the *PCAtools* (Blighe and Lun 2020) package for the purpose
of performing PCA.

``` r
  library(PCAtools)
  p <- pca(assay(sce, 'scaled'), metadata = metadata(sce))

  biplot(p,
    lab = NULL,
    pointSize = 1.0,
    colby = 'treatment',
    legendPosition = 'right',
    title = 'PCA applied to CyTOF data',
    caption = paste0('25000 cells randomly selected after ',
      'having filtered for low variance'))
```

![Perform PCA](README_files/figure-gfm/ex1-1.png)

We can add the rotated component loadings as a new reduced dimensional
component to our dataset. Let’s just add the first 20 PCs.

``` r
  reducedDim(sce, 'PCA') <- p$rotated[,1:20]
```

For more functionality via *PCAtools*, check the vignette: [PCAtools:
everything Principal Component
Analysis](https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

## Perform UMAP

UMAP can be performed on the entire dataset, if your computer’s memory
will permit. Currently it’s default is to use the data contained in the
‘scaled’ assay component of your *SingleCellExperiment* object.

``` r
  sce <- performUMAP(sce)
```

UMAP can also be stratified based on a column in your metadata, e.g.,
(treated versus untreated samples); however, to do this, I recommend
creating separate *SingleCellExperiment* objects from the very start,
i.e., from the the data input stage, and processing the data separately
for each group.

We can also perform UMAP on a select number of PC eigenvectors.
*PCAtools* (Blighe and Lun 2020) can be used to infer ideal number of
dimensions to use via the elbow method and Horn’s parallel analysis.

``` r
  elbow <- findElbowPoint(p$variance)
  horn <- parallelPCA(assay(sce, 'scaled'))

  elbow
```

    ## PC9 
    ##   9

``` r
  horn$n
```

    ## [1] 4

Let’s use the number of PCs identified by Horn’s.

``` r
  sce <- performUMAP(sce, reducedDim = 'PCA', dims = c(1:horn$n))
```

## Create a contour plot of the UMAP layout

This and the remaining sections in this tutorial are about producing
great visualisations of the data and attempting to make sense of it,
while not fully overlapping with functionalioty provided by other
programs that operate in tis space.

With the contour plot, we are essentially looking at celluar density. It
can provide for a beautiful viusualisation in a manuscript while also
serving as a useful QC tool: if the density is ‘scrunched up’ into a
single area in the plot space, then there are likely issues with your
input data distribution. We want to see well-separated, high density
‘islands’, or, at least, gradual gradients that blend into one another
across high density ‘peaks’.

``` r
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
```

![Create a contour plot of the UMAP
layout](README_files/figure-gfm/ex2-1.png)

## Show marker expression across the layout

Here, we randomly select some markers and then plot their expression
profiles across the UMAP layouts.

``` r
  markers <- sample(rownames(sce), 6)
  markers
```

    ## [1] "CD25" "CD2"  "CD58" "CD22" "CD19" "CD47"

``` r
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
```

![Show marker expression across the
layout](README_files/figure-gfm/ex3-1.png)

## Shade cells by metadata

Shading cells by metadata can be useful for identifying any batch
effects, but also useful for visualising, e.g., differences across
treatments.

First, let’s take a look inside the metadata that we have.

``` r
  head(metadata(sce))
```

    ##       sample   group treatment
    ## cell1    P04 Disease       CD3
    ## cell2  HD262 Healthy    Unstim
    ## cell3  HD262 Healthy    Unstim
    ## cell4    P00 Disease    Unstim
    ## cell5  HD262 Healthy      CD46
    ## cell6  HD262 Healthy    Unstim

``` r
  levels(metadata(sce)$group)
```

    ## [1] "Healthy" "Disease"

``` r
  levels(metadata(sce)$treatment)
```

    ## [1] "CD46"   "Unstim" "CD3"

``` r
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
```

![Shade cells by metadata](README_files/figure-gfm/ex4-1.png)

## Find ideal clusters in the UMAP layout via k-nearest neighbours

This function utilises the k nearest neighbours (k-NN) approach from
Seurat, which works quite well on flow cytometry and CyTOF UMAP layouts,
from my experience.

``` r
  sce <- clusKNN(sce,
    k.param = 20,
    prune.SNN = 1/15,
    resolution = 0.01,
    algorithm = 2,
    verbose = FALSE)

  sce <- clusKNN(sce,
    reducedDim = 'UMAP_PCA',
    clusterAssignName = 'Cluster_PCA',
    k.param = 20,
    prune.SNN = 1/15,
    resolution = 0.01,
    algorithm = 2,
    verbose = FALSE)

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
```

![Find ideal clusters in the UMAP layout via k-nearest
neighbours](README_files/figure-gfm/ex5-1.png)

## Plot marker expression per identified cluster

``` r
  markerExpressionPerCluster(sce,
    caption = 'Cluster assignments based on UMAP performed on expression values',
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)
```

![Plot marker expression per identified
cluster1](README_files/figure-gfm/ex6a-1.png)

``` r
  clusters <- unique(metadata(sce)[['Cluster_PCA']])
  clusters
```

    ## [1] 1 2 6 0 5 3 4 7

``` r
  markers <- sample(rownames(sce), 8)
  markers
```

    ## [1] "CD42" "CD62" "CD8"  "CD1"  "CD20" "CD53" "CD36" "CD34"

``` r
  markerExpressionPerCluster(sce,
    clusters = clusters,
    clusterAssign = metadata(sce)[['Cluster_PCA']],
    markers = markers,
    nrow = 2, ncol = 5,
    caption = 'Cluster assignments based on UMAP performed on PC eigenvectors',
    stripLabSize = 22,
    axisLabSize = 22,
    titleLabSize = 22,
    subtitleLabSize = 18,
    captionLabSize = 18)
```

![Plot marker expression per identified
cluster2](README_files/figure-gfm/ex6b-1.png)

Try all markers across a single cluster:

``` r
  cluster <- sample(unique(metadata(sce)[['Cluster']]), 1)
  cluster
```

    ## [1] 7

``` r
  markerExpressionPerCluster(sce,
    clusters = cluster,
    markers = rownames(sce),
    stripLabSize = 20,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 14,
    captionLabSize = 12)
```

![Plot marker expression per identified
cluster3](README_files/figure-gfm/ex6c-1.png)

## Determine enriched markers in each cluster and plot the expression signature

This method also calculates metacluster abundances across a chosen
phenotype. The function returns a data-frame, which can then be exported
to do other analyses.

### Disease vs Healthy metacluster abundances

``` r
  markerEnrichment(sce,
    method = 'quantile',
    studyvarID = 'group')
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:right;">

Cluster

</th>

<th style="text-align:right;">

nCells

</th>

<th style="text-align:right;">

TotalCells

</th>

<th style="text-align:right;">

PercentCells

</th>

<th style="text-align:left;">

NegMarkers

</th>

<th style="text-align:left;">

PosMarkers

</th>

<th style="text-align:right;">

PerCent\_HD00

</th>

<th style="text-align:right;">

PerCent\_HD01

</th>

<th style="text-align:right;">

PerCent\_HD262

</th>

<th style="text-align:right;">

PerCent\_P00

</th>

<th style="text-align:right;">

PerCent\_P02

</th>

<th style="text-align:right;">

PerCent\_P03

</th>

<th style="text-align:right;">

PerCent\_P04

</th>

<th style="text-align:right;">

PerCent\_P08

</th>

<th style="text-align:right;">

nCell\_Healthy

</th>

<th style="text-align:right;">

nCell\_Disease

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

5902

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

23.608

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD21+CD23+CD30+CD32+CD37+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

63.3852931

</td>

<td style="text-align:right;">

0.1016605

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

36.5130464

</td>

<td style="text-align:right;">

3741

</td>

<td style="text-align:right;">

2161

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

4451

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

17.804

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD9-CD10-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

22.6915300

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0224669

</td>

<td style="text-align:right;">

0.0449337

</td>

<td style="text-align:right;">

18.5800944

</td>

<td style="text-align:right;">

11.817569

</td>

<td style="text-align:right;">

46.8434060

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

1011

</td>

<td style="text-align:right;">

3440

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3877

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

15.508

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

0.0257931

</td>

<td style="text-align:right;">

32.9120454

</td>

<td style="text-align:right;">

0.0773794

</td>

<td style="text-align:right;">

66.8042301

</td>

<td style="text-align:right;">

0.0257931

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.1289657

</td>

<td style="text-align:right;">

0.0257931

</td>

<td style="text-align:right;">

1280

</td>

<td style="text-align:right;">

2597

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2972

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

11.888

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD9-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD32+CD39+CD47+CD57+CD58+

</td>

<td style="text-align:right;">

0.0672948

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

20.9623149

</td>

<td style="text-align:right;">

30.551817

</td>

<td style="text-align:right;">

48.4185734

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2970

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

2900

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

11.600

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD38+CD47+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0344828

</td>

<td style="text-align:right;">

77.9655172

</td>

<td style="text-align:right;">

0.3448276

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.0689655

</td>

<td style="text-align:right;">

21.5862069

</td>

<td style="text-align:right;">

2262

</td>

<td style="text-align:right;">

638

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

2028

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

8.112

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

45.7593688

</td>

<td style="text-align:right;">

0.0493097

</td>

<td style="text-align:right;">

53.9940828

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.1479290

</td>

<td style="text-align:right;">

0.0493097

</td>

<td style="text-align:right;">

929

</td>

<td style="text-align:right;">

1099

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

1895

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

7.580

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD11-CD15-CD55-

</td>

<td style="text-align:left;">

CD1+CD21+CD23+CD31+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

1.1609499

</td>

<td style="text-align:right;">

0.2110818

</td>

<td style="text-align:right;">

98.3641161

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.2110818

</td>

<td style="text-align:right;">

0.0527704

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

1869

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

915

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

3.660

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD30+CD32+CD37+CD49+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

58.1420765

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

41.8579235

</td>

<td style="text-align:right;">

532

</td>

<td style="text-align:right;">

383

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

0.148

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD12-CD17-CD44-CD55-

</td>

<td style="text-align:left;">

CD1+CD32+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

100.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

0.092

</td>

<td style="text-align:left;">

CD2-CD4-CD6-CD7-CD8-CD9-CD10-CD11-CD18-CD55-

</td>

<td style="text-align:left;">

CD1+CD36+CD41+CD47+CD57+CD58+

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

13.0434783

</td>

<td style="text-align:right;">

4.347826

</td>

<td style="text-align:right;">

82.6086957

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

23

</td>

</tr>

</tbody>

</table>

.

### Treatment type metacluster abundances

``` r
  markerEnrichment(sce,
    sampleAbundances = FALSE,
    method = 'quantile',
    studyvarID = 'treatment')
```

v

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:right;">

Cluster

</th>

<th style="text-align:right;">

nCells

</th>

<th style="text-align:right;">

TotalCells

</th>

<th style="text-align:right;">

PercentCells

</th>

<th style="text-align:left;">

NegMarkers

</th>

<th style="text-align:left;">

PosMarkers

</th>

<th style="text-align:right;">

nCell\_CD46

</th>

<th style="text-align:right;">

nCell\_Unstim

</th>

<th style="text-align:right;">

nCell\_CD3

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

5902

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

23.608

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD21+CD23+CD30+CD32+CD37+CD57+CD58+

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

5900

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

4451

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

17.804

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD9-CD10-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

1012

</td>

<td style="text-align:right;">

3421

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3877

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

15.508

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

3875

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2972

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

11.888

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD9-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD32+CD39+CD47+CD57+CD58+

</td>

<td style="text-align:right;">

2923

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

49

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

2900

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

11.600

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD38+CD47+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

2900

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

2028

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

8.112

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

19

</td>

<td style="text-align:right;">

2006

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

1895

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

7.580

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD8-CD11-CD15-CD55-

</td>

<td style="text-align:left;">

CD1+CD21+CD23+CD31+CD32+CD57+CD58+

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

1883

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

915

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

3.660

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD15-

</td>

<td style="text-align:left;">

CD1+CD23+CD30+CD32+CD37+CD49+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

912

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

0.148

</td>

<td style="text-align:left;">

CD2-CD4-CD7-CD11-CD12-CD17-CD44-CD55-

</td>

<td style="text-align:left;">

CD1+CD32+CD51+CD57+CD58+

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

25000

</td>

<td style="text-align:right;">

0.092

</td>

<td style="text-align:left;">

CD2-CD4-CD6-CD7-CD8-CD9-CD10-CD11-CD18-CD55-

</td>

<td style="text-align:left;">

CD1+CD36+CD41+CD47+CD57+CD58+

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

11

</td>

</tr>

</tbody>

</table>

.

### Expression signature

The expression signature is a quick way to visualise which markers are
more or less expressed in each identified cluster of cells.

``` r
  plotSignatures(sce,
    labCex = 1.5,
    legendCex = 1.5)
```

![Determine enriched markers in each cluster and plot the expression
signature](README_files/figure-gfm/ex7-1.png)

# Tutorial 2: Import from Seurat

Due to the fact that *scDataviz* is based on *SingleCellExperiment*, it
has increased interoperability with other packages, including the
popular *Seurat* (Stuart et al. 2018). Taking the data produced from the
[Seurat
Tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) on
Peripheral Blood Mononuclear Cells (PBMCs), we can convert this to a
*SingleCellExperiment* object recognisable by *scDataviz* via
`as.SingleCellExperiment()`.

For a full workflow for *Seurat*-to-*scDataviz*, see the [GitHub
vignette](https://github.com/kevinblighe/scDataviz/).

# Tutorial 3: Import any numerical data

*scDataviz* will work with any numerical data, too. Here, we show a
quick example of how one can import a data-matrix of randomly-generated
numbers that follow a negative binomial distribution, comprising 2500
cells and 20 markers:

``` r
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
```

    ##       group
    ## cell1     A
    ## cell2     A
    ## cell3     A
    ## cell4     A
    ## cell5     A
    ## cell6     A

``` r
  sce <- importData(mat,
    assayname = 'normcounts',
    metadata = metadata)
  sce
```

    ## class: SingleCellExperiment 
    ## dim: 20 2500 
    ## metadata(1): group
    ## assays(1): normcounts
    ## rownames(20): CD1 CD2 ... CD19 CD20
    ## rowData names(0):
    ## colnames(2500): cell1 cell2 ... cell2499 cell2500
    ## colData names(0):
    ## reducedDimNames(0):
    ## spikeNames(0):
    ## altExpNames(0):

This will also work without any assigned metadata.

``` r
  sce <- importData(mat,
    assayname = 'normcounts',
    metadata = NULL)
  sce
```

    ## class: SingleCellExperiment 
    ## dim: 20 2500 
    ## metadata(0):
    ## assays(1): normcounts
    ## rownames(20): CD1 CD2 ... CD19 CD20
    ## rowData names(0):
    ## colnames(2500): cell1 cell2 ... cell2499 cell2500
    ## colData names(0):
    ## reducedDimNames(0):
    ## spikeNames(0):
    ## altExpNames(0):

# Acknowledgments

  - Shahram Kordasti
  - Jessica Timms
  - James Opzoomer
  - Marcel Ramos (Bioconductor)
  - Bioinformatics CRO

# Session info

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/atlas-base/atlas/libblas.so.3.0
    ## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=pt_BR.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=pt_BR.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] PCAtools_2.1.4              cowplot_1.0.0              
    ##  [3] lattice_0.20-41             reshape2_1.4.4             
    ##  [5] ggrepel_0.8.2               ggplot2_3.3.0              
    ##  [7] scDataviz_0.99.48           SingleCellExperiment_1.8.0 
    ##  [9] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
    ## [11] BiocParallel_1.20.1         matrixStats_0.56.0         
    ## [13] Biobase_2.46.0              GenomicRanges_1.38.0       
    ## [15] GenomeInfoDb_1.22.1         IRanges_2.20.2             
    ## [17] S4Vectors_0.24.4            BiocGenerics_0.32.0        
    ## [19] kableExtra_1.1.0            knitr_1.28                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] corrplot_0.84            sn_1.6-1                 plyr_1.8.6              
    ##   [4] igraph_1.2.5             lazyeval_0.2.2           splines_3.6.3           
    ##   [7] flowCore_1.52.1          listenv_0.8.0            TH.data_1.0-10          
    ##  [10] digest_0.6.25            htmltools_0.4.0          gdata_2.18.0            
    ##  [13] fansi_0.4.1              magrittr_1.5             cluster_2.1.0           
    ##  [16] ROCR_1.0-7               globals_0.12.5           readr_1.3.1             
    ##  [19] sandwich_2.5-1           askpass_1.1              colorspace_1.4-1        
    ##  [22] rvest_0.3.5              rappdirs_0.3.1           xfun_0.13               
    ##  [25] dplyr_0.8.5              crayon_1.3.4             RCurl_1.98-1.2          
    ##  [28] jsonlite_1.6.1           survival_3.1-12          zoo_1.8-7               
    ##  [31] ape_5.3                  glue_1.4.0               gtable_0.3.0            
    ##  [34] zlibbioc_1.32.0          XVector_0.26.0           webshot_0.5.2           
    ##  [37] leiden_0.3.3             BiocSingular_1.2.2       future.apply_1.4.0      
    ##  [40] scales_1.1.0             mvtnorm_1.1-0            bibtex_0.4.2.2          
    ##  [43] Rcpp_1.0.4.6             isoband_0.2.1            metap_1.3               
    ##  [46] plotrix_3.7-7            viridisLite_0.3.0        dqrng_0.2.1             
    ##  [49] reticulate_1.15          rsvd_1.0.3               tsne_0.1-3              
    ##  [52] umap_0.2.5.0             htmlwidgets_1.5.1        httr_1.4.1              
    ##  [55] gplots_3.0.3             RColorBrewer_1.1-2       TFisher_0.2.0           
    ##  [58] ellipsis_0.3.0           Seurat_3.1.4             ica_1.0-2               
    ##  [61] farver_2.0.3             pkgconfig_2.0.3          uwot_0.1.8              
    ##  [64] labeling_0.3             tidyselect_1.0.0         rlang_0.4.5             
    ##  [67] munsell_0.5.0            tools_3.6.3              cli_2.0.2               
    ##  [70] ggridges_0.5.2           evaluate_0.14            stringr_1.4.0           
    ##  [73] yaml_2.2.1               npsurv_0.4-0             fitdistrplus_1.0-14     
    ##  [76] caTools_1.18.0           purrr_0.3.3              RANN_2.6.1              
    ##  [79] pbapply_1.4-2            future_1.16.0            nlme_3.1-147            
    ##  [82] xml2_1.3.1               compiler_3.6.3           rstudioapi_0.11         
    ##  [85] plotly_4.9.2.1           png_0.1-7                lsei_1.2-0              
    ##  [88] tibble_3.0.0             stringi_1.4.6            highr_0.8               
    ##  [91] RSpectra_0.16-0          Matrix_1.2-18            multtest_2.42.0         
    ##  [94] vctrs_0.2.4              mutoss_0.1-12            pillar_1.4.3            
    ##  [97] lifecycle_0.2.0          Rdpack_0.11-1            lmtest_0.9-37           
    ## [100] RcppAnnoy_0.0.16         data.table_1.12.8        bitops_1.0-6            
    ## [103] irlba_2.3.3              gbRd_0.4-11              patchwork_1.0.0         
    ## [106] R6_2.4.1                 KernSmooth_2.23-16       gridExtra_2.3           
    ## [109] codetools_0.2-16         MASS_7.3-51.5            gtools_3.8.2            
    ## [112] assertthat_0.2.1         openssl_1.4.1            withr_2.1.2             
    ## [115] sctransform_0.2.1        mnormt_1.5-6             multcomp_1.4-13         
    ## [118] GenomeInfoDbData_1.2.2   hms_0.5.3                grid_3.6.3              
    ## [121] tidyr_1.0.2              DelayedMatrixStats_1.8.0 rmarkdown_2.1           
    ## [124] Rtsne_0.15               numDeriv_2016.8-1.1

# References

﻿Blighe (2020)

Blighe and Lun (2020)

Lun and Risso (2020)

Stuart et al. (2018)

Wickham (2016)

<div id="refs" class="references">

<div id="ref-scDataviz">

Blighe, K. 2020. “scDataviz: single cell dataviz and downstream
analyses.” <https://github.com/kevinblighe/scDataviz.>

</div>

<div id="ref-PCAtools">

Blighe, K, and A Lun. 2020. “PCAtools: everything Principal Component
Analysis.” <https://github.com/kevinblighe/PCAtools.>

</div>

<div id="ref-Lun">

Lun, A, and D Risso. 2020. “SingleCellExperiment: S4 Classes for Single
Cell Data.” <https://bioconductor.org/packages/SingleCellExperiment.>

</div>

<div id="ref-satijalab">

Stuart, Tim, Andrew Butler, Paul Hoffman, Christoph Hafemeister,
Efthymia Papalexi, William M Mauck III, Marlon Stoeckius, Peter Smibert,
and Rahul Satija. 2018. “Comprehensive Integration of Single Cell Data.”
*bioRxiv*. <https://doi.org/10.1101/460147>.

</div>

<div id="ref-Wickham">

Wickham, H. 2016. “ggplot2: Elegant Graphics for Data Analysis.”
Springer-Verlag New York, ISBN: 978-3-319-24277-4.

</div>

</div>
