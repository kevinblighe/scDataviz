contourplot <- function(
  data)
{
  gg <- ggplot(as.data.frame(data$layout), aes(UMAP1, UMAP2)) +

    stat_density2d(aes(alpha = ..level.., fill = ..level..), size = 1, bins = 350, geom = "polygon") +

    scale_fill_gradient(low = "darkblue", high = "darkred") +

    scale_alpha(range = c(0.00, 0.5), guide = FALSE) +

    geom_density2d(colour="black") +

    #geom_point(data = mat, size = 0.01, shape = ".") +

    #Set the size of the plotting window
    theme_bw(base_size=24) +

    #Modify various aspects of the plot text and legend
    theme(
      legend.position = "right",
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 26, face = "bold", vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 26, face = "bold", vjust = 1),
      plot.caption = element_text(angle = 0, size = 26, face = 'plain', vjust = 1),

      axis.text.x = element_text(angle = 0, size = 26, face = "bold", hjust = 0.5),
      axis.text.y = element_text(angle = 0, size = 26, face = "bold", vjust = 0.5),
      axis.title = element_text(size = 26, face = "bold"),

      #Legend
      legend.key = element_blank(),     #removes the border
      legend.key.size = unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text = element_text(size = 26),  #Text size
      title = element_text(size = 26)) +      #Title text size

    #Change the size of the icons/symbols in the legend
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    #guides(alpha=FALSE)

    #Set x- and y-axes labels
    xlab("UMAP 1") +
    ylab("UMAP 2") +

    #xlim(-15, 10) +
    #ylim(-13, 14) +

  labs(
    title = "",
    subtitle = "Cellular density and contours",
    caption = paste0("Total cells, ", nrow(data), "; Bins, 350"))

  return(gg)
}