basetheme <- function(
  titleLabSize,
  subtitleLabSize,
  captionLabSize,
  axisLabSize,
  xlabAngle,
  xlabhjust,
  xlabvjust,
  ylabAngle,
  ylabhjust,
  ylabvjust,
  legendPosition,
  legendLabSize) {

  theme_bw(base_size = 24) +

  theme(
    legend.background = element_rect(),

    title = element_text(size = legendLabSize),

    plot.title = element_text(angle = 0, size = titleLabSize,
      face = 'bold', vjust=1),
    plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
      face = 'plain', vjust = 1),
    plot.caption = element_text(angle = 0, size = captionLabSize,
      face = 'plain', vjust = 1),

    axis.text.x = element_text(angle = xlabAngle, size = axisLabSize,
      hjust = xlabhjust, vjust = xlabvjust),
    axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
      hjust = ylabhjust, vjust = ylabvjust),
    axis.title = element_text(size = axisLabSize),

    legend.title = element_blank(),
    legend.position = legendPosition,
    legend.key = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = legendLabSize))
}
