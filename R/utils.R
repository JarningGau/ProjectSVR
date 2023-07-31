#' @importFrom magrittr "%<>%"
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' multi-core FindAllMarkers()
#' @param seu Seurat object
#' @param do.flatten Flattern the marker list into data.frame. Default: TRUE
#' @param only.pos Only positive markers. Default: TRUE
#' @param n.cores The number of threads. Default: 10
#' @param ... Arguments passed to Seurat::FindMarkers()
#' @return data.frame or list contains putative markers.
#' @references https://gist.github.com/diazdc/1735102c243cd16acb1b1f3fd09a26e1
#' @concept extract_signature
#' @export
mcFindAllMarkers <- function(seu, do.flatten=TRUE, only.pos=TRUE, n.cores=10, ...) {
  n_clust <- levels(Seurat::Idents(seu))
  find.markers <- function(i){
    ident1 <- i
    ident2 <- n_clust[n_clust != i]
    table <- Seurat::FindMarkers(seu, ident.1 = ident1, ident.2 = ident2, only.pos = only.pos, ...)
    table$Gene.name.uniq <- rownames(table)
    table$cluster <- rep(i, nrow(table))
    return(table)
  }

  marker_results <- parallel::mclapply(n_clust, find.markers, mc.cores = n.cores)
  names(marker_results) <- n_clust
  # nice way to flatten list into a single DF
  if (do.flatten) {
    marker_results <- dplyr::bind_rows(marker_results)
  }
  return(marker_results)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Internal ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Run sklearn.neighbors.NearestNeighbors
#'
#' @param data Input data
#' @param query Data to query against data
#' @param k Number of nearest neighbors to compute
#' @param ... other parameters pass to \code{sklearn.neighbors.NearestNeighbors()}.
#' @concept internal
#'
sknn <- function(data, query=data, k=5L, ...) {
  args <- list(...)
  args <- args[intersect(names(args), c("radius", "algorithm", "leaf_size", "metric", "p", "n_jobs"))]
  if ("leaf_size" %in% names(args)) {
    args$leaf_size %<>% as.integer()
  }
  if ("p" %in% names(args)) {
    args$p %<>% as.integer()
  }
  if ("n_jobs" %in% names(args)) {
    args$n_jobs %<>% as.integer()
  }
  args$n_neighbors <- as.integer(k)
  nn <- reticulate::import("sklearn.neighbors")
  knn <- do.call(what = nn$NearestNeighbors, args = args)
  knn <- knn$fit(as.matrix(data))
  results <- knn$kneighbors(as.matrix(query))
  names(results) <- c("nn.dists", "nn.idx")
  # Note: python is 0-indices while R is 1-indices.
  results$nn.idx <- results$nn.idx + 1
  return(results)
}


#' Internal helper function to dispatch to various neighbor finding methods
#'
#' @param data Input data
#' @param query Data to query against data
#' @param k Number of nearest neighbors to compute
#' @param method "rann", "sklearn"
#' @param ... Other args pass \code{nn2()} or \code{sknn()}.
#' @concept internal
#' @importFrom RANN nn2
NNHelper <- function(data, query = data, k, method, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = nn2, args = args)
      },
      "sklearn" = {
        ## This codes is not neccessary because sknn will check the parameters.
        # args <- args[intersect(x = names(x = args), y = names(x = formals(fun = sknn)))]
        do.call(what = 'sknn', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'sklearn'")
    )
  )
  n.ob <- list(
    nn.idx = results$nn.idx, # (n, k)
    nn.dists = results$nn.dists, # (n, k)
    ref.cell.names = rownames(x = data),
    query.cell.names = rownames(x = query) # length = n
  )
  class(n.ob) <- "Neighbor"
  return(n.ob)
}


# Query the nearest cells for each query cell
#
# @param object Neighbor object.
#
# @return The k nearest cell names for each query cell.
cells.Neighbor <- function(object) {
  ref.cell.names <- object$ref.cell.names
  nn.idx <- object$nn.idx # (n, k) (cells, knn)
  # return nn cells except itself
  results <- lapply(1:nrow(nn.idx), function(xx) ref.cell.names[nn.idx[xx, -1]] )
  names(results) <- object$query.cell.names
  return(results)
}


cosine.dist <- function(X) {
  sim <- X / sqrt(rowSums(X * X))
  sim <- sim %*% t(sim)
  D_sim <- 1 - sim
  return(D_sim)
}


# L1 normalization
L1norm <- function(A) {
  A <- as.matrix(A)
  C <- matrix(rep(rowSums(A), ncol(A)), ncol = ncol(A))
  C <- A / C
  colnames(C) <- colnames(A)
  C
}


# L2 Normalization
L2norm <- function(A) {
  A <- as.matrix(A)
  C <- matrix(rep(sqrt(rowSums(A * A)), ncol(A)), ncol = ncol(A))
  C <- A / C
  colnames(C) <- colnames(A)
  C
}

