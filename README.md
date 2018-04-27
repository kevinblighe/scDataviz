# cytofNet
Unbiased identification of ideal clusters, k, in mass cytometry time-of-flight (CyTOF) data, with interrogation of clusters via network analysis.
<h2>Tutorial</h2>
This tutorial covers just 3 samples (sample1.fcs, sample2.fcs, sample3.fcs) that have already had the basic QC performed on them, e.g., removal of dead cells and manual gating (although, manual gating is supported).
Key points:

 - data is downsampled based on low variance

 - unbiased clustering to identify ideal number of centers, k, is performed via bootstrapped partitioning around medoids (PAM) and 3 metrics: Gap statistic; silhouette coefficient; elbow method

 - medoids from PAM are ued to infer low/high marker expression

 - network plot construction used to show relationships between identified clusters

<h3>1, Setup / initialisation</h3>

```{r}
  #Set CPU cores for parallel-related functions
  cpucores <- 16
  require(parallel)
  options("mc.cores"=cpucores)

  #Set CPU cores for doParallel-related functions
  require(doParallel)
  cores <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cores)
```

<h3>2, Set global variables</h3>

```{r}
  #Set background noise threshold - values below this are set to 0
  BackgroundNoiseThreshold <- 1

  #Euclidean norm threshold - this is the square root of the sum of all the squares
  EuclideanNormThreshold <- 1

  #Choose a transformation function (any mathematical function)
  transFun <- function (x) asinh(x)

  #Set hyperbolic arc-sine factor (NB - asinh(x/5) is recommended for CyTOF and FACS data)
  asinhFactor <- 5
```

<h3>3, Data input and conversion (FCS -> CSV)</h3>

```{r}
  #Convert FCS to CSV
  require(flowCore)
  source("R/fcs2csv.R")
  fcs2csv(sample1.fcs, sample1.csv)
  fcs2csv(sample2.fcs, sample2.csv)
  fcs2csv(sample3.fcs, sample3.csv)
  
  sample1 <- read.csv("sample1.csv")
  sample2 <- read.csv("sample2.csv")
  sample3 <- read.csv("sample3.csv")

  #Create a vector of all variable names for the samples
  AllSamples <- c("sample1", "sample2", "sample3")
```

<h3>4, Histograms to check distribution of data</h3>

```{r}
  par(mfrow=c(1,3), cex=1.2)
  source("R/transform.R")
  colours <- c("red", "gold", "skyblue")
  for (i in 1:length(AllSamples))
  {
    x <- as.matrix(get(AllSamples[1]))
    x <- x[,-which(colnames(x) %in% c("DNA.1", "DeadLive"))]
    x <- transform(x, BackgroundNoiseThreshold, EuclideanNormThreshold, transFun, asinhFactor)
    hist(data.matrix(x), main=paste("Hyperbolic arc-sine\nsample",  i), breaks=30, col=colours[i])
  }
```
<img src="images/checkDistribution.png"></img>

<h3>5, Hexagonal binning to capture and summarise the 3-dimensional nature of CyTOF data</h3>

```{r}
require(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(9,"BuPu")))
source("R/facsplot.R")
facsplot("CD8", "HLA.DR", get(AllSamples[1]), bins=400, main="FACS plot", xlim=c(-10,10), ylim=c(-10,10), x1=0, x2=5, y1=0, y2=5, cex=1.0, colramp=rf) # HexBin ignores NAs and doesn't count them whilst binning
```
<img src="images/facsplot.png"></img>

<hr>

<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (University College London)</li>
  <li>Kevin Blighe (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Davide Lucchesi (Queen Mary University of London)</li>
</ul>
