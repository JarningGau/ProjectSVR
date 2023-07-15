# Merge cells into metacells to speed up cNMF.
#' @include utils.R
#' @importFrom methods setClass new
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Class definitions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The GridDensity class
#'
#' The GridDensity class is an intermediate data class that storing the density of
#' cells in each bin (also each mesh points) and the cell assignments of the nearest
#' mesh point for performing the downstream analysis - production of pseudo mini bulk
#' RNA-seq counts and grid differential expression analysis.
#'
#' @slot data a matrix records embedding matrix.
#' @slot mesh.points a table records coordinates (primary key, <X,Y>), knn.density,
#' bin coordinates, bin.ID, and bin.density of each mesh point in data.frame.
#' @slot cell2mp a named vector, records cell assignments to the nearest mesh point.
#' @concept class
#'
GridDensity <- setClass(
  Class = "GridDensity",
  slots = list(
    data = "matrix",
    mesh.points = "data.frame",
    cell2mp = "vector"
  )
)


#' subset a GridDensity object.
#'
#' @param x, GridDensity object
#' @param bin.density.threshold Only bins with density greater than this value are retained. Default: mean-3*sd
#' @param mesh.density.threshold Only mesh points with density greater than this value are retained.
#' We do not filter on mesh point directly. Default: -1
#' @param n.cells.threshold Only mesh points with assigned cells no less than this value are retained. Default: 10
#' @param n.mesh.threshold Only bins with number of mesh pints no less than this value are retained. Default: -1
#' @param ... Arguments passed to other methods.
#' Default: max number of mesh points per bin minus 1.
#' @concept meta_cell
#' @export
subset.GridDensity <- function(x, bin.density.threshold=NULL, mesh.density.threshold=-1, n.cells.threshold=10, n.mesh.threshold=-1, ...) {
  mesh.info <- x@mesh.points
  bin.info <- mesh.info %>%
    dplyr::select(.data$bin.ID, .data$bin.density) %>%
    dplyr::distinct()
  if (is.null(bin.density.threshold)) {
    bin.density.threshold <- mean(bin.info$bin.ID) - 3*stats::sd(bin.info$bin.ID)
  }
  if (is.null(n.mesh.threshold)) {
    n.mesh.threshold <- max(table(mesh.info$bin.ID)) - 1
  }
  n.cells <- table(x@cell2mp)
  # 1. filter on bins: bin.density > threshold
  # 2. filter on mesh points: mesh.density > threshold
  mesh.info.new <- mesh.info %>%
    dplyr::filter(.data$bin.density > bin.density.threshold) %>%
    dplyr::filter(.data$mesh.density > mesh.density.threshold)

  # 3. filter on mesh points: number of assigned cells >= threshold
  kept.mesh <- names(n.cells[n.cells >= n.cells.threshold])
  mesh.info.new <- mesh.info.new %>% dplyr::filter(.data$mesh.ID %in% kept.mesh)

  # 4. filter on bins: number of mesh points >= threshold
  n.mesh <- table(mesh.info.new$bin.ID)
  kept.bins <- names(n.mesh[n.mesh >= n.mesh.threshold])
  mesh.info.new <- dplyr::filter(mesh.info.new, .data$bin.ID %in% kept.bins)
  # 5. filter cell2mp
  cell2mp <- x@cell2mp
  cell2mp <- cell2mp[cell2mp %in% mesh.info.new$mesh.ID]
  x@mesh.points <- mesh.info.new
  x@cell2mp <- cell2mp
  return(x)
}


#' Grid density plot
#'
#' @param x A \code{GridDensity} object
#' @param ... Arguments passed to other methods.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_tile scale_fill_viridis_c labs theme_bw theme element_blank element_text
#' @importFrom ggplot2 geom_density
#' @concept meta_cell
plot.GridDensity <- function(x, ...) {
  object <- x
  cell.emb <- data.frame(X = object@data[,1], Y = object@data[,2])
  bin.data <- object@mesh.points %>%
    dplyr::select(.data$bin.X, .data$bin.Y, .data$bin.density) %>%
    dplyr::distinct()
  p1 <- ggplot() +
    geom_point(aes(.data$X, .data$Y), data=cell.emb, size=0.1, alpha=0.2, color="grey", show.legend=F) +
    # geom_raster(aes(.data$bin.X, .data$bin.Y, fill=.data$bin.density), data=bin.data) +
    geom_tile(aes(.data$bin.X, .data$bin.Y, fill=.data$bin.density), data=bin.data) +
    scale_fill_viridis_c() +
    labs(x = colnames(object@data)[1], y = colnames(object@data)[2]) +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank())
  p2 <- ggplot(bin.data) +
    geom_density(aes(.data$bin.density), fill="red", color="red", alpha=0.1) +
    labs(x = "Bin density") +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  cowplot::plot_grid(p1,p2, nrow = 1, rel_widths = c(7,4))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Internal functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Sum the densities of each point over the bins
#'
#' @param x matrix, n.col=2
#' @param bins matrix, n.col=2
#' @param mesh.density length(mesh.density) equals nrow(x). Default: 1
#' @return data.frame of the mesh points. (X,Y) position, bin.ID, mesh.density, bin.density.
#' @keywords internal
.hist2d <- function(x, bins, mesh.density=1) {
  data.tmp <- data.frame(
    mesh.ID = rownames(x),
    mesh.X = x[,1],
    mesh.Y = x[,2],
    mesh.density = mesh.density,
    bin.ind.X = findInterval(x[,1], bins[,1]),
    bin.ind.Y = findInterval(x[,2], bins[,2])
  ) %>%
    dplyr::mutate(bin.X = bins[.data$bin.ind.X, 1], bin.Y = bins[.data$bin.ind.Y, 2]) %>%
    dplyr::mutate(bin.ID = paste(.data$bin.ind.X, .data$bin.ind.Y, sep = "_"), .before = "bin.X") %>%
    dplyr::select(-.data$bin.ind.X, -.data$bin.ind.Y)

  # sum the whole grid should be 1
  bin.density <- data.tmp %>%
    dplyr::select(.data$bin.ID, .data$mesh.density) %>%
    dplyr::group_by(.data$bin.ID) %>%
    dplyr::summarise(bin.density = mean(.data$mesh.density))
  bin.density <- bin.density %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin.density = .data$bin.density / sum(.data$bin.density))

  dplyr::left_join(data.tmp, bin.density, by="bin.ID")
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### API ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Estimate KNN density on 2d embedding.
#' @importFrom magrittr "%<>%"
#'
#' @param emb.mat a embedding matrix, rows are cells, columns are dimensions.
#' @param k Number of neighbors. Default: 10
#' @param n.bins number of bins for density resampling. Default: 50
#' @param n.mesh in each bin, density will be calculated around (n.mesh^2) points. Default: 3
#' @param do.plot whether plot results. Default: TRUE.
#' @returns a GridDensity object
#' @concept meta_cell
#' @export
EstimateKnnDensity <- function(emb.mat, k=10, n.bins=50, n.mesh=2, do.plot=TRUE) {
  # check parameters
  if (ncol(emb.mat) != 2) {
    stop("Error dimensions of embedding matrix (emb.mat). Only support 2d coordinates.")
  }
  if (is.double(k) || is.double(n.bins) || is.double(n.mesh)) {
    k %<>% as.integer()
    n.bins %<>% as.integer()
    n.mesh %<>% as.integer()
  }
  X <- emb.mat[, 1]
  Y <- emb.mat[, 2]
  # Ref to knnDREMI() - https://github1s.com/KrishnaswamyLab/scprep/blob/master/scprep/stats.py
  # 1. create bin and mesh points
  x.bins <- seq(min(X), max(X), length.out = n.bins)
  y.bins <- seq(min(Y), max(Y), length.out = n.bins)
  x.mesh <- seq(min(X), max(X), length.out = n.mesh*n.bins)
  y.mesh <- seq(min(Y), max(Y), length.out = n.mesh*n.bins)

  # 2. create mesh points
  mesh.points <- cbind(
    rep(x.mesh, each = length(y.mesh)),
    rep(y.mesh, times = length(x.mesh))
  )
  x.mesh.ind <- findInterval(mesh.points[,1], x.mesh)
  y.mesh.ind <- findInterval(mesh.points[,2], y.mesh)
  rownames(mesh.points) <- paste(x.mesh.ind, y.mesh.ind, sep="_")
  colnames(mesh.points) <- colnames(emb.mat)

  # 3. create a kNN graph between cells and mesh points.
  nn.ranked <- NNHelper(data = emb.mat, query = mesh.points, k = k, method = "rann")

  # 4. calculate the kNN density around the mesh points.
  # Get area, density of each point
  area = pi * (nn.ranked$nn.dists[, k] ^ 2)
  density = k / area
  # get list of all mesh points that are not bin intersections
  # Note: the difference between `|` and `||` is `|` is vecterized boolean operator.
  mesh.mask <- (mesh.points[, 1] %in% x.bins) | (mesh.points[, 2] %in% y.bins)
  mesh.info <- .hist2d(x = mesh.points[!mesh.mask, ], bins = cbind(x.bins, y.bins), mesh.density = density[!mesh.mask])

  # 5. assign cells to nearest mesh.points
  nn.ranked <- NNHelper(data = mesh.points[!mesh.mask, ], query = emb.mat, k = 1, method = "rann")
  indices <- as.vector(nn.ranked$nn.idx)

  cell2mp <- rownames(mesh.points[!mesh.mask, ])[indices]
  names(cell2mp) <- rownames(emb.mat)
  object <- new("GridDensity", data = emb.mat,  mesh.points = mesh.info, cell2mp = cell2mp)
  if (do.plot) {
    print(plot(object))
  }
  return(object)
}


#' Merge cells around mesh points in each bin into pseudo-bulk gene expression matrix.
#'
#' @param dge Matrix or sparse matrix for gene expression matrix rows are genes, columns are cells.
#' @param gd A GridDensity object.
#' @param by Merge raw counts by "bins" or "mesh.points". Default: "mesh.points"
#' @param cores Number of threads, -1 means all available threads. Default: -1
#' @return A gene expression matrix, rows are genes, columns are merged cells.
#' @concept meta_cell
#' @export
MergeCells <- function(dge, gd, by = "mesh.points", cores=-1) {
  # check parameters
  cell2mp <- gd@cell2mp
  outliers <- setdiff(names(cell2mp), colnames(dge))
  if (length(outliers) > 0) {
    stop(sprintf("Outlier cells in GridDensity object: %s\n Please check your inputs.",
         paste(outliers, collapse = ",")))
  }
  if (!by %in% c("bins", "mesh.points")) {
    stop("Invalid 'by'. Please choose one of 'bins', 'mesh.points'")
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())

  mesh.points <- gd@mesh.points
  mesh2bin <- mesh.points$bin.ID
  names(mesh2bin) <- mesh.points$mesh.ID
  cellmeta <- data.frame(
    cell.ID = names(cell2mp),
    mesh.ID = cell2mp,
    bin.ID = mesh2bin[cell2mp]
  )
  bin.IDs <- unique(cellmeta$bin.ID)
  dge.bin.list <- parallel::mclapply(bin.IDs, function(xx) {
    cells.use <- dplyr::filter(cellmeta, .data$bin.ID == xx)$cell.ID
    dge[, cells.use]
  }, mc.cores = cores)
  names(dge.bin.list) <- bin.IDs

  if (by == "mesh.points") {
    message("Merging cells by mesh points ...")
    mesh.IDs <- unique(cellmeta$mesh.ID)
    bulk.dge <- pbapply::pbsapply(mesh.IDs, function(xx) {
      bin.use <- mesh2bin[xx]
      cells.use <- dplyr::filter(cellmeta, .data$mesh.ID == xx)$cell.ID
      Matrix::rowSums(dge.bin.list[[bin.use]][, cells.use])
    })
  } else if (by == "bins") {
    bulk.dge <- pbapply::pbsapply(dge.bin.list, Matrix::rowSums)
  } else {
    stop("Invalid 'by'. Please choose one of 'bins', 'mesh.points'")
  }
  return(bulk.dge)
}


