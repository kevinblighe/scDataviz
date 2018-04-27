# cytofNet
Unbiased identification of ideal clusters, k, in mass cytometry time-of-flight (CyTOF) data, with interrogation of clusters via network analysis.
<h1>Tutorial</h1>
<code>
  #Remove scientific notation and set decimal places to 3
  options(scipen=100, digits=3)

  #Set CPU cores for parallel-related functions
  cpucores <- 16
  require(parallel)
  options("mc.cores"=cpucores)

  #Set CPU cores for doParallel-related functions
  require(doParallel)
  cores <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cores)

  #Convert FCS to CSV
  require(flowCore)
  source("testing/fcs2csv.R")
  fcs2csv(sample1.fcs, sample1.csv)
  fcs2csv(sample2.fcs, sample2.csv)
  fcs2csv(sample3.fcs, sample3.csv)
  
  sample1 <- read.csv("sample1.csv")
  sample2 <- read.csv("sample2.csv")
  sample3 <- read.csv("sample3.csv")

  #Create a vector of all variable names for the samples
  AllSamples <- c("sample1", "sample2", "sample3")
</code>
<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (University College London)</li>
  <li>Kevin Blighe (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Davide Lucchesi (Queen Mary University of London)</li>
</ul>
