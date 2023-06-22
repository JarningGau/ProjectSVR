#' @include themes.R
NULL


#' Visualizing the projected query cells on reference via a density plot.
#'
#' @importFrom ggplot2 aes geom_point geom_text scale_color_manual guides facet_wrap
#' @importFrom rlang .data
#'
#' @param seu.q A query Seurat object containing the projected reference embeddings.
#' @param reference Reference model.
#' @param ref.color.by Name of column in reference data to be used for coloring points.
#' @param ref.colors Vector of colors to be used for reference data points. Default: NULL
#' @param ref.alpha Alpha value for reference data points. Default: 0.5
#' @param ref.size Size of reference data points. Default: 0.5
#' @param ref.emb.col The column names storing the reference embeddings in reference model.
#' @param query.reduction Reduction name to use for query data. Default: 'ref.umap'
#' @param query.size Size of query data points. Default: 1
#' @param query.alpha Alpha value for query data points. Default: 1
#' @param split.by Column name to split query data by. Default: NULL
#' @param n.row Number of rows for facet_wrap when splitting by column. Default: NULL
#' @param legend.ncol Number of columns for legend. Default: 1
#'
#' @return A ggplot object.
#' @export
#'
PlotProjection <- function(seu.q, reference, ref.color.by,
                           ref.colors = NULL, ref.alpha = 0.5, ref.size = 0.5,
                           ref.emb.col = paste0("UMAP_", 1:2),
                           query.reduction = "ref.umap",
                           query.size = 1, query.alpha = 1, split.by = NULL,
                           n.row = NULL, legend.ncol = 1) {
  # Check if seu.q is a Seurat object
  if (!inherits(seu.q, "Seurat")) {
    stop("seu.q argument must be a Seurat object")
  }
  # Check if reference is a named list with ref.cellmeta and ref.assay
  if (!is.list(reference) || !all(c("ref.cellmeta") %in% names(reference))) {
    stop("reference argument must be a named list with ref.cellmeta element")
  }
  # Check if ref.color.by is a valid column name in reference$ref.cellmeta$meta.data
  if (!is.null(ref.color.by) && !ref.color.by %in% colnames(reference$ref.cellmeta$meta.data)) {
    stop(paste0("Invalid ref.color.by argument '", ref.color.by, "'. It must be a valid column name in reference$ref.cellmeta$meta.data"))
  }
  # Check if ref.colors is NULL or a valid color vector of length equal to the number of levels in ref.color.by
  if (!is.null(ref.colors) && (!is.vector(ref.colors) || length(ref.colors) != nlevels(reference$ref.cellmeta$meta.data[, ref.color.by]))) {
    stop("Invalid ref.colors argument. It must be NULL or a valid color vector of length equal to the number of levels in ref.color.by")
  }
  #
  if (!all(ref.emb.col %in% colnames(reference$ref.cellmeta$meta.data)) || length(ref.emb.col) != 2) {
    stop(paste0("Invalid ref.emb.col argument '", paste(ref.emb.col), "'. They must be valid column names in reference$ref.cellmeta$meta.data"))
  }
  # Check if query.reduction is a valid reduction name in seu.q
  if (!query.reduction %in% names(seu.q@reductions)) {
    stop(paste0("Invalid query.reduction argument '", query.reduction, "'. It must be a valid reduction name in seu.q"))
  }
  # Check if split.by is NULL or a valid column name in query.data
  if (!is.null(split.by) && !split.by %in% colnames(seu.q@meta.data)) {
    stop(paste0("Invalid split.by argument '", split.by, "'. It must be a valid column name in query.data"))
  }

  ref.data <- reference$ref.cellmeta$meta.data
  ref.text <- reference$ref.cellmeta$text.pos
  if (is.null(ref.colors)) {
    ref.colors <- reference$ref.cellmeta$colors
  }
  ref.p <- ggplot(ref.data, aes(get(ref.emb.col[1]), get(ref.emb.col[2]), color = get(ref.color.by))) +
    geom_point(size = ref.size, alpha = ref.alpha) +
    xlab("Dimension 1") + ylab("Dimension 2")

  if (!is.null(ref.text)) {
    ref.p <- ref.p +
      geom_text(inherit.aes = F, data = ref.text, mapping = aes(.data$x, .data$y, label = .data$label), size = 5)
  }
  if (!is.null(ref.colors)) {
    ref.p <- ref.p + scale_color_manual(values = ref.colors)
  }
  ref.p <- ref.p +
    guides(color = guide_legend(title = "", override.aes = list(size = 3, alpha = 1),
                                ncol = legend.ncol)) +
    DimTheme()

  dim.key <- seu.q[[query.reduction]]@key
  qx <- paste0(dim.key, 1)
  qy <- paste0(dim.key, 2)
  query.data <- Seurat::FetchData(seu.q, vars = c(qx, qy, split.by))

  p <- ref.p +
    geom_point(inherit.aes = F, data = query.data,
               mapping = aes(get(qx), get(qy)), shape = 17, alpha = query.alpha, size = query.size) +
    geom_density_2d(inherit.aes = F, data = query.data,
                    mapping = aes(get(qx), get(qy)), color = "blue", contour_var="ndensity")
  if (!is.null(split.by)) {
    if (is.null(n.row)) {
      n.row <- floor(sqrt(length(unique(query.data[[split.by]])))) + 1
    }
    p <- p + facet_wrap(~get(split.by), nrow = n.row)
  }
  p
}

#' Calculate cell population percentage
#'
#' This function calculates percentage statistics based on a given dataset and grouping variables.
#'
#' @param cellmeta A data frame containing cell metadata.
#' @param by A character vector specifying the variable(s) used to group the data.
#' @param fill A character vector specifying the variable used to fill the data.
#'
#' @return A data frame containing the calculated percentage statistics.
#'
#' @examples
#' PercentageStat(mtcars, "cyl", "gear")
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
PercentageStat <- function(cellmeta, by, fill) {
  # Check if cellmeta is a data.frame or tibble
  if (!is.data.frame(cellmeta)) {
    stop("cellmeta must be a data.frame or tibble.")
  }

  # Check if by and fill are character vectors
  if (!is.character(by) || !is.character(fill)) {
    stop("by and fill must be character vectors.")
  }

  # Check if by and fill exist in cellmeta
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta))) {
    stop("by and fill must be columns in cellmeta.")
  }
  data.stat <- as.data.frame(table(cellmeta[[by]], cellmeta[[fill]]))
  colnames(data.stat)[1:2] <- c(by, fill)
  data.stat <- data.stat %>%
    group_by_at(by) %>%
    mutate(margin.freq = sum(.data$Freq)) %>%
    mutate(proportion = .data$Freq / .data$margin.freq)
  return(data.stat)
}

#' Plot Alluvial Plot
#'
#' @param cellmeta a data frame containing the cell metadata.
#' @param by a character string specifying the column name in \code{cellmeta} which will be used to group the cells.
#' @param fill a character string specifying the column name in \code{cellmeta} which will be used to fill the plot.
#' @param colors (optional) vector of colours to use for filling the plot. If not specified, the default colour scheme will be used.
#' @param bar.width a numeric value between 0 and 1 specifying the width of the bars in the plot.
#' @param legend.ncol an integer specifying the number of columns in the legend.
#'
#' @return A ggplot object representing an alluvial plot.
#'
#' @importFrom ggplot2 ggplot geom_area geom_col scale_x_continuous scale_y_continuous guides guide_legend
#' @importFrom dplyr group_by arrange mutate ungroup summarise filter
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' AlluviaPlot(mtcars, "cyl", "gear")
#'
#' @export
#' @references https://stackoverflow.com/questions/73372641/shaded-area-between-bars-using-ggplot2
AlluviaPlot <- function(cellmeta, by, fill, colors = NULL, bar.width = 0.5, legend.ncol = 1) {
  pop.stat <- PercentageStat(cellmeta, by, fill)

  alluvia <- pop.stat %>%
    group_by(get(by)) %>%
    arrange(get(by), get(fill), by_group = TRUE) %>%
    mutate(y = .data$proportion) %>%
    ungroup() %>%
    mutate(x = as.numeric(as.factor(get(by)))) %>%
    group_by(get(by), get(fill)) %>%
    summarise(x = c(.data$x - 0.25, .data$x, .data$x + 0.25), y = .data$y)
  colnames(alluvia)[1:2] <- c(by, fill)

  x.labels <- unique(alluvia[[by]])

  p <- ggplot(alluvia %>% filter(.data$x %% 1 == 0), aes(.data$x, .data$y, fill = get(fill))) +
    geom_area(data = alluvia, alpha = 0.4, position = 'fill') +
    geom_col(width = bar.width, color = 'gray50', position = 'fill') +
    scale_x_continuous(breaks = 1:length(x.labels), labels = x.labels) +
    scale_y_continuous(labels = scales::label_percent()) +
    guides(fill = guide_legend(ncol = legend.ncol))

  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  p <- p +
    theme_classic(base_size = 15) +
    theme(legend.title = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  p
}

#' Create Group Preference Plot
#'
#' A function to create a heatmap that displays the O/E ratio of a two-column table grouped by a selected column.
#'
#' @param cellmeta A data frame.
#' @param group.by A character string indicating the name of the column in `cellmeta` used for grouping.
#' @param preference.on A character string indicating the name of the column in `cellmeta` used for calculating O/E ratio.
#' @param palette A color palette to use in the plot. Default is "Blues".
#' @param column_names_rot degree (0-360) to rotate x-axis labels. Default is 0.
#' @param ... additional arguments passed to the 'Heatmap' function from ComplexHeatmap package.
#'
#' @return A heatmap that displays the O/E ratio of a two-column table grouped by a selected column.
#'
#' @import ggplot2
#' @import ComplexHeatmap
#' @importFrom grid grid.text gpar
#' @importFrom stats quantile
#'
#' @examples
#' data(mtcars)
#' mtcars$cyl <- factor(mtcars$cyl)
#' mtcars$carb <- factor(mtcars$carb)
#' GroupPreferencePlot(mtcars, "cyl", "carb")
#'
#' @export
GroupPreferencePlot <- function(cellmeta, group.by, preference.on, palette = "Blues", column_names_rot = 0, ...) {
  # Check if cellmeta is a data frame
  if (!is.data.frame(cellmeta)) {
    stop("cellmeta must be a data frame")
  }

  # Check if group.by and preference.on are character variables in cellmeta
  if (!(group.by %in% names(cellmeta)) || !(preference.on %in% names(cellmeta))) {
    stop("group.by and preference.on must be character variables in cellmeta")
  }

  ## calculate O/E ratio
  pop.stat.mat <- table(cellmeta[[group.by]], cellmeta[[preference.on]])
  pop.stat.mat.o <- apply(pop.stat.mat, 1, function(xx) xx / sum(xx)) %>% t()
  pop.stat.mat.e <- colSums(pop.stat.mat) / sum(pop.stat.mat)
  pop.stat.mat.e <- matrix(rep(pop.stat.mat.e, nrow(pop.stat.mat)), ncol = nrow(pop.stat.mat)) %>% t()
  pop.stat.mat.oe <- t(pop.stat.mat.o / pop.stat.mat.e)
  pop.stat.mat.oe[is.na(pop.stat.mat.oe)] <- 1
  col.names.order <- levels(cellmeta[[group.by]])
  row.names.order <- rev(levels(cellmeta[[preference.on]]))
  if (is.null(row.names.order)) {
    pop.stat.mat.oe <- pop.stat.mat.oe[order(pop.stat.mat.oe[,1], decreasing = T), ]
  } else {
    pop.stat.mat.oe <- pop.stat.mat.oe[row.names.order, ]
  }
  if (is.null(col.names.order)) {
    pop.stat.mat.oe <- pop.stat.mat.oe[, order(pop.stat.mat.oe[1,], decreasing = T)]
  } else {
    pop.stat.mat.oe <- pop.stat.mat.oe[, col.names.order]
  }

  ## plots

  cut.off <- quantile(pop.stat.mat.oe, .95, na.rm = TRUE)
  cut.off2 <- round(quantile(pop.stat.mat.oe, .9, na.rm = TRUE), 1)
  pop.stat.mat.oe.tile <- ifelse(pop.stat.mat.oe > cut.off, cut.off, pop.stat.mat.oe)

  colors <- tryCatch({
    RColorBrewer::brewer.pal(n = 7, name = palette)
  }, error = function(e) {
    message(e$message, "Use default: Blues")
    RColorBrewer::brewer.pal(n = 7, name = "Blues")
  })

  Heatmap(pop.stat.mat.oe.tile, name = "Ratio (o/e)", cluster_rows = F, cluster_columns = F, col = colors,
          column_names_rot = column_names_rot,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if (pop.stat.mat.oe[i, j] > cut.off2) {
              grid.text(sprintf("%.2f", pop.stat.mat.oe[i, j]), x, y, gp = gpar(fontsize=10, col="white", fontface="bold"))
            } else {
              grid.text(sprintf("%.2f", pop.stat.mat.oe[i, j]), x, y, gp = gpar(fontsize=10, col="black"))
            }
          }, ...)
}

#' Differential Cellular Abundance Test
#'
#' Computes a Wilcoxon rank sum test for each cell type between two groups.
#'
#' @param cellmeta A data frame containing metadata information about cells.
#' @param celltype.col The name of the column in \code{cellmeta} that contains cell type information.
#' @param sample.col The name of the column in \code{cellmeta} that contains sample identifiers.
#' @param group.col The name of the column in \code{cellmeta} that specifies the grouping variable.
#'
#' @return A data frame containing results of Wilcoxon rank sum tests for each cell type.
#'
#' @details This function takes in a data frame \code{cellmeta}, a character string \code{celltype.col},
#' a character string \code{sample.col}, and a character string \code{group.col} as input. It calculates
#' fold change, mean percentage, p-value for each cell type using Wilcoxon rank sum test.
#'
#' @examples
#' library(dplyr)
#' data(iris)
#'
#' # Add some metadata to iris data frame
#' iris$Sample_ID <- rep(paste0("sample", 1:10), each = 15)
#' iris$Group <- rep(c("A","B"), each = 75)
#'
#' # Run the AbundanceTest function
#' AbundanceTest(cellmeta = iris,
#' celltype.col = "Species",
#' sample.col = "Sample_ID",
#' group.col = "Group")
#'
#' @importFrom stats wilcox.test
#' @importFrom magrittr %>%
#' @importFrom dplyr distinct
#'
#' @export
AbundanceTest <- function(cellmeta, celltype.col, sample.col, group.col) {
  stopifnot(is.data.frame(cellmeta),
            is.character(celltype.col),
            is.character(sample.col),
            is.character(group.col))

  # Check that celltype.col exists in cellmeta
  if (!(celltype.col %in% colnames(cellmeta))) {
    stop(paste0(celltype.col, " not found in cellmeta"))
  }

  # Check that sample.col exists in cellmeta
  if (!(sample.col %in% colnames(cellmeta))) {
    stop(paste0(sample.col, " not found in cellmeta"))
  }

  # Check that group.col exists in cellmeta
  if (!(group.col %in% colnames(cellmeta))) {
    stop(paste0(group.col, " not found in cellmeta"))
  }
  count.matrix <- table(cellmeta[[celltype.col]], cellmeta[[sample.col]])
  count.matrix <- as.matrix(count.matrix)
  count.matrix.norm <- apply(count.matrix, 2, function(xx) xx / sum(xx))

  sample.meta <- cellmeta[, c(sample.col, group.col)] %>%
    distinct() %>%
    as.data.frame()
  rownames(sample.meta) <- sample.meta[[sample.col]]
  sample.meta[[sample.col]] <- NULL
  sample.meta <- sample.meta[colnames(count.matrix), ]

  celltypes <- sort(rownames(count.matrix.norm))

  group.levels <- levels(cellmeta[[group.col]])
  if (is.null(group.levels)) {
    group.levels <- unique(cellmeta[[group.col]])
  }

  test.df <- lapply(celltypes, function(ct) {
    c1 <- count.matrix.norm[ct, sample.meta == group.levels[1]]
    c2 <- count.matrix.norm[ct, sample.meta == group.levels[2]]
    res <- wilcox.test(c1, c2)
    data.frame(
      celltype = ct,
      fold.change = mean(c2) / mean(c1),
      p.value = res$p.value,
      mean.perc = mean(c(c1,c2))
    )
  })
  test.df <- do.call(rbind, test.df)
  return(test.df)
}

#' Plot the Volcano Plot
#'
#' Volcano plot for visualizing the cellular abundance difference. This function
#' accepts the outputs of `AbundanceTest()`.
#'
#' @param test.df A data frame from `AbundanceTest()`.
#' @param p.value.cutoff The p-value cutoff for showing the cell type labels.
#' Default: 0.001
#' @param xlab The label for the x-axis. Default: NULL
#' @param ylab The label for the y-axis. Default: NULL
#' @param colors A vector of colors to be used in the plot.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export
VolcanoPlot <- function(test.df, p.value.cutoff = 0.001, xlab = NULL, ylab = NULL, colors = NULL) {
  p <- ggplot(test.df, aes(log2(.data$fold.change), .data$mean.perc,
                           color = .data$celltype, size = -log10(.data$p.value) )) +
    geom_point() +
    geom_vline(xintercept = c(-1,1), color = "blue", linetype = "dashed") +
    ggrepel::geom_text_repel(
      inherit.aes = F, data = subset(test.df, .data$p.value < p.value.cutoff),
      mapping = aes(log2(.data$fold.change), .data$mean.perc, label = .data$celltype),
      nudge_x = .1, nudge_y = .01, size = 4) +
    scale_y_continuous(labels = scales::percent_format(1L)) +
    theme_classic(base_size = 15) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "top")
  if (!is.null(xlab)) {
    p <- p + xlab(label = xlab)
  }
  if (!is.null(ylab)) {
    p <- p + ylab(label = ylab)
  }
  if (is.null(colors)) {
    # use the default color pallate
    colors <- scales::hue_pal()(length(unique(test.df$celltype)))
  }
  p <- p + scale_color_manual(values = colors, guide= "none")
  p
}
