#' Package-wide, non-user function used to set a base \code{ggplot2} theme.
#'
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param xlabAngle Rotation angle of x-axis labels.
#' @param xlabhjust Horizontal adjustment of x-axis labels.
#' @param xlabvjust Vertical adjustment of x-axis labels.
#' @param ylabAngle Rotation angle of y-axis labels.
#' @param ylabhjust Horizontal adjustment of y-axis labels.
#' @param ylabvjust Vertical adjustment of y-axis labels.
#' @param legendPosition Position of \code{legend ('top', 'bottom',
#'   'left', 'right', 'none')}.
#' @param legendLabSize Size of plot legend text.
#'
#' @details
#' Package-wide, non-user function used to set a base \code{ggplot2} theme.
#'
#' @return A \code{list} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#' # create a theme
#' th <- basetheme(
#'   titleLabSize = 16,
#'   subtitleLabSize = 12,
#'   captionLabSize = 12,
#'   axisLabSize = 16,
#'   xlabAngle = 0,
#'   xlabhjust = 0.5,
#'   xlabvjust = 0.5,
#'   ylabAngle = 0,
#'   ylabhjust = 0.5,
#'   ylabvjust = 0.5,
#'   legendPosition = 'none',
#'   legendLabSize = 12)
#'
#' @importFrom ggplot2 theme_bw theme
#'
#' @export
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
