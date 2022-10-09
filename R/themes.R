#' @import ggplot2
NULL

#' theme for DimPlot
#' @export
DimTheme <- function() {
  theme_bw(base_size = 15) %+replace%
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = .5, face = "bold"),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}
