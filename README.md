# cytofNet
Unbiased identification of ideal clusters, k, in mass cytometry time-of-flight (CyTOF) data, with interrogation of clusters via network analysis.
<h2>Tutorial</h2>

<h3>Setup / initialisation</h3>

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

<h3>Set global variables</h3>

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

<h3>Data input and conversion (FCS -> CSV)</h3>

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

<h3>Histograms to check distribution of data</h3>

```{r}
  par(mfrow=c(1,3), cex=1.2)
  source("R/transform.R")
  x <- as.matrix(get(AllSamples[1]))
  x <- x[,-which(colnames(x) %in% c("DNA.1", "DeadLive"))]
  x <- transform(x, BackgroundNoiseThreshold, EuclideanNormThreshold, transFun, asinhFactor)
  hist(data.matrix(x), main="Hyperbolic arc-sine\nsample 1", breaks=30, col="red")
  hist(data.matrix(x), main="Hyperbolic arc-sine\nsample 2", breaks=30, col="gold")
  hist(data.matrix(x), main="Hyperbolic arc-sine\nsample 3", breaks=30, col="skyblue")
  <img src="images/checkDistribution.png"></img>
```

<hr>

<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (University College London)</li>
  <li>Kevin Blighe (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Davide Lucchesi (Queen Mary University of London)</li>
</ul>
