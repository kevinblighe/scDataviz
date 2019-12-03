medianHeatmap <- function(
  data,
  feature = NULL
)
{
  require(pheatmap)
  require(RColorBrewer)
  require(dplyr)
  
  # Get the median marker expression per sample
  data.median = data.frame(sample_id = data@metadata$group, t(assay(sct))) %>%
    group_by(sample_id) %>% 
    summarize_all(list(median))
  
  # reshape to matrix
  data.median.sample = t(data.median[, -1])
  colnames(data.median.sample) = data.median$sample_id
  
  # Using metadata feature slot
  if(is.null(feature) == FALSE){
    
    # Summarize sampleID to metadata feature 
    feature.summary = data.frame(sample_id = data@metadata$group, feature = data@metadata[feature]) %>%
      distinct()
    
    # Match feature
    mm = match(colnames(data.median.sample), feature.summary$sample_id)
    annotation_col = data.frame(feature = feature.summary[, feature][mm],
                               row.names = colnames(data.median.sample))
    
    # possible colors
    anno_colors =  brewer.pal(n = 9, name = "Set1")
     
    # Create annotation 
    condition = anno_colors[1:length(unique(feature.summary[,feature]))]
    names(condition) = unique(feature.summary[,feature])                      
    annotation_colors = list(feature = condition)
    
    # Colors for the heatmap
    color_bar = brewer.pal(n = 9, name = "Greens")
    
    # execute heatmap
    heatmap = pheatmap(data.median.sample, color = color_bar, 
                        display_numbers = F,
                        fontsize_number = 9, 
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors, 
                        clustering_method = "average")
    
    return(heatmap)
    
  } else{
    
    # Colors for the heatmap
    color_bar = brewer.pal(n = 9, name = "Greens")
    
    # execute heatmap
    heatmap = pheatmap(data.median.sample, 
                       color = color_bar, 
                       display_numbers = F,
                      clustering_method = "average")
    
    return(heatmap)
    
  }
  
  
  
  
}
