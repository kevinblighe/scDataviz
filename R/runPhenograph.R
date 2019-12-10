runPhenograph <- function(
  sce,
  k = 30,
  markers = NULL
){
  # Need to implement marker selection
  require(Rphenograph)
  
  if(is.null(feature) == FALSE){
    
    data = t(assay(sce))
    
    rphenograph_output = Rphenograph(data[,markers], k = k)
  
    identity = factor(membership(rphenograph_output[[2]]))
  
    sce@metadata$phenograph = identity
  
  return(sce)
  
  }
  else{
    
    rphenograph_output <- Rphenograph(t(assay(sce)))
    
    identity <- factor(membership(rphenograph_output[[2]]))
    
    sce@metadata$phenograph = identity
    
    return(sce)
    
  }
  
}