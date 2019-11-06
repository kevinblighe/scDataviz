plotNRS <- function(
  data,
  feature = NULL
  )
  {
  require(reshape2)
  require(ggplot2)
  
  ## Define a function that calculates the NRS per sample
  ## Function taken from Levine et al.
  NRS = function(x, ncomp = 3){
    pr =prcomp(x, center = TRUE, scale. = FALSE)
    score =rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
    return(score)
  }
  
  # Split the data into a list by sample
  data.list = data.frame(sample_id = data$metadata$group, data$expression) %>% 
              group_by(sample_id) %>% group_split(keep = FALSE) 
  
  # Summarize sampleID to metadata feature 
  feature.summary = data.frame(sample_id = data$metadata$group, feature = data$metadata[feature]) %>%
    distinct()
  
  # calculate the score
  nrs_sample = sapply(data.list, NRS)
  colnames(nrs_sample) = feature.summary$sample_id
  nrs = rowMeans(nrs_sample, na.rm = TRUE)
  
  ## Plot the NRS for ordered markers
  nrs_ord_markers = names(sort(nrs, decreasing = TRUE))
  nrs_sample = data.frame(t(nrs_sample))
  nrs_sample$sample_id = feature.summary$sample_id
  
  # melt for plotting
  nrs_boxplot = melt(nrs_sample, id.var = "sample_id",
               value.name = "nrs", variable.name = "antigen")
  
  nrs_boxplot$antigen = factor(nrs_boxplot$antigen, levels = nrs_ord_markers)
  
  # Using metadata feature slot
  if(is.null(feature) == FALSE){
    
    # match the condition
    mm = match(nrs_boxplot$sample_id, feature.summary$sample_id)
    nrs_boxplot$condition = feature.summary$condition[mm]
  
    ncol = length(unique(nrs_boxplot$antigen))
  
    # Plot the boxplot
    gbox = ggplot(nrs_boxplot, aes(x = condition, y = nrs)) +
            geom_point(aes_string(color = feature), alpha = 0.9,
                      position = position_jitter(width = 0.5, height = 0)) +
            geom_boxplot(aes(color = condition), outlier.color = NA, fill = NA, lwd=1.5) +
            theme_classic() +
            facet_wrap('antigen', scales = 'fix', ncol = ncol) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    
    return(gbox)
  
  } else{
    
    # Plot the boxplot
    gbox = ggplot(nrs_boxplot, aes(x = antigen, y = nrs)) +
      geom_point(alpha = 0.9, position = position_jitter(width = 0.3, height = 0)) +
      geom_boxplot(outlier.color = NA, fill = NA) +
      stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    
    return(gbox)
      
    
  }
  
  
}
