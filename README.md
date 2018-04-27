# cytofNet
Unbiased identification of ideal clusters, k, in mass cytometry time-of-flight (CyTOF) data, with interrogation of clusters via network analysis.
<h1>Tutorial</h1>
<code>#Remove scientific notation and set decimal places to 3</code>
<code>options(scipen=100, digits=3)</code>

<code>#Set CPU cores for parallel-related functions</code>
<code>cpucores <- 16</code>
<code>require(parallel)</code>
<code>options("mc.cores"=cpucores)</code>

<code>#Set CPU cores for doParallel-related functions</code>
<code>require(doParallel)</code>
<code>cores <- makeCluster(detectCores(), type='PSOCK')</code>
<code>registerDoParallel(cores)</code>

<code>#Convert FCS to CSV</code>
<code>require(flowCore)</code>
<code>source("testing/fcs2csv.R")</code>
<code>fcs2csv(sample1.fcs, sample1.csv)</code>
<code>fcs2csv(sample2.fcs, sample2.csv)</code>
<code>fcs2csv(sample3.fcs, sample3.csv)</code>
  
<code>sample1 <- read.csv("sample1.csv")</code>
<code>sample2 <- read.csv("sample2.csv")</code>
<code>sample3 <- read.csv("sample3.csv")</code>

<code>#Create a vector of all variable names for the samples</code>
<code>AllSamples <- c("sample1", "sample2", "sample3")</code>
<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (University College London)</li>
  <li>Kevin Blighe (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Kevin Blighe (Queen Mary University of London)</li>
  <li>Davide Lucchesi (Queen Mary University of London)</li>
</ul>
