setClusterIdents <- function(
  sce,
  orig.ident,
  new.ident,
  clustering = 'phenograph'
){
  
  ## Convert to factor with merged clusters in desired order
  # This would mean I have to order all the factors, there are too many
  levels_clusters_merged = unique(new.ident)
  new.ident = factor(new.ident, levels = levels_clusters_merged) 
  
  # Map the clusters to a new annotation
  mm <- match(sce@metadata[clustering][,1], orig.ident)
  cell_annotation = new.ident[mm]
  
  sce@metadata$cell_annotation = cell_annotation
  
  return(sce)
  
}