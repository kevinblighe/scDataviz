mdsplot <- function(
  data,
  feature = NULL
  )
  {
  require(dplyr)
  require(ggplot2)
  require(limma)
  require(ggrepel)
  
  # Get the median marker expression per sample 
  data.median = data.frame(sample_id = data$metadata$group, data$expression) %>%
                group_by(sample_id) %>% 
                summarize_all(list(median))
  
  data.median.sample = t(data.median[, -1])
  
  colnames(data.median.sample) = data.median$sample_id
  mds = plotMDS(data.median.sample, plot = FALSE)
  mdsdf = data.frame(MDS1 = mds$x, MDS2 = mds$y,
                    sample_id = colnames(data.median.sample))
  
  # Using metadata feature slot
  if(is.null(feature) == FALSE){
    
    # Summarize sampleID to metadata feature 
    feature.summary = data.frame(sample_id = data$metadata$group, feature = data$metadata[feature]) %>%
                      distinct()
    
    # Match feature
    mm = match(mdsdf$sample_id, feature.summary$sample_id)
    mdsdf[,feature] = as.factor(feature.summary[,feature][mm])
  
    # Plot
    mds = ggplot(mdsdf, aes_string(x = "MDS1", y = "MDS2", color = feature)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_label_repel(aes(label = sample_id)) +
    theme_classic() +
    coord_fixed()
    
    return(mds)
  
  } else{
    
    # Plot
    mds = ggplot(mdsdf, aes(x = MDS1, y = MDS2)) +
        geom_point(size = 2, alpha = 0.8) +
        geom_label_repel(aes(label = sample_id)) +
        theme_classic() +
        coord_fixed()
    
    return(mds)
    
  }
  
  }
  