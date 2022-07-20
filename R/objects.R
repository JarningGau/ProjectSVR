#' The CellProject class
#'
#' The main class used by Monocle to hold single cell expression data. CellDataSet extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from a matrix of expression values Methods that operate on CellDataSet objects constitute the basic Monocle workflow.
#'
#' @slot embeddings Matrix containing the predicted embeddings computed by ProjectNewdata.
#' @slot refined.embeddings Matrix containing the refined predicted embeddings computed by RefineProjection
#' @slot data data.frame containing the source data for cell projection.
#' @slot cellmeta data.frame containing the data generated during the cell projection analysis.
#' @slot neighbors \code{Neighbor} class containing the nearest neighbors for calculating the projection metrics.
#'
#' @name CellProject-class
#' @rdname CellProject-class
#' @aliases CellProject-class
#' @exportClass CellProject
setClass(
  Class = "CellProject",
  slots = c(
    embeddings = "matrix",
    refined.embeddings = "matrix",
    data = "data.frame",
    cellmeta = "data.frame",
    neighbors = "ANY")
)


#' Creates a new CellProject object.
#'
#' @param embeddings matrix containing predicted embeddings
#' @param data matrix containing source data for cell projection
#' @param cellmeta data frame containing cell meta data generated during cell projection
#' @return a new CellProject object
newCellProject <- function(embeddings, data, cellmeta = NULL){
  if (missing(embeddings) || missing(data)) {
    stop("Must provide 'embeddings' and 'data'")
  }
  if (is.null(rownames(data))) {
    stop("No cell names (rownames) present in the input data")
  }
  if (is.null(colnames(data))) {
    stop("No features names (colnames) present in the input data")
  }
  if (is.null(cellmeta)) {
    cellmeta <- data.frame(row.names = rownames(data))
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  if (!is.matrix(embeddings)) {
    embeddings <- as.matrix(embeddings)
  }
  if (!all(rownames(embeddings) == rownames(data))) {
    stop("The rownames of 'embeddings' and 'data' are not matched!")
  }
  object <- new(
    Class = "CellProject",
    embeddings = embeddings,
    refined.embeddings = new(Class = "matrix"),
    data = data,
    cellmeta = cellmeta,
    neighbors = NULL)
  # validObject(object)
  return(object)
}


# setValidity("CellProject", function(object) {
#   emb <- object@embeddings
#   data <- object@data
#   if (!all(rownames(emb) == rownames(data))) {
#     "The embeddings and source data are not matched!"
#   }
#   TRUE
# } )


#' Merge CellProject Objects
#'
#' @param x object
#' @param y another object
#' @param ... arguments to be passed to or from methods.
#' @return merged object
#' @method merge CellProject
#'
#' @export
merge.CellProject <- function(x, y, ...) {
  embeddings <- rbind(x@embeddings, y@embeddings)
  data <- rbind(x@data, y@data)
  cellmeta <- rbind(x@cellmeta, y@cellmeta)

  check.dup.cn <- duplicated(colnames(cellmeta))
  if (any(check.dup.cn)) {
    warning("Duplicated colnames in 'cellmeta' will auto dropped.")
    cellmeta <- cellmeta[!check.dup.cn]
  }
  merged.object <- newCellProject(
    embeddings = embeddings,
    data = data,
    cellmeta = cellmeta)
  if (nrow(x@refined.embeddings) > 0 &&
      nrow(y@refined.embeddings) > 0 &&
      ncol(x@refined.embeddings) == ncol(y@refined.embeddings)) {
    refined.embeddings <- rbind(x@refined.embeddings, y@refined.embeddings)
    merged.object@refined.embeddings <- refined.embeddings
  }
  return(merged.object)
}


#' Subset CellProject Object
#' @param x A CellProject object
#' @param cells cells for subsetting
#'
#' @return A subsetted CellProject object
#'
subset2 <- function(x, cells) {
  if (missing(x) || missing(cells)) {
    stop("Must provide 'x' and 'cells'")
  }
  object <- newCellProject(
    embeddings = x@embeddings[cells, ],
    data = x@data[cells, ],
    cellmeta = x@cellmeta[cells, , drop=FALSE]
  )
  if (nrow(x@refined.embeddings) > 0) {
    object@refined.embeddings <- x@refined.embeddings[cells, ]
  }
  return(object)
}


#' Split CellProject Object
#' @param x CellProject object
#' @param split.by Attribute for splitting
#'
#' @return A list of CellProject object
#'
#' @export
SplitCellProject <- function(x, split.by) {
  if (missing(x) || missing(split.by)) {
    stop("Must provide 'x' and 'split.by'")
  }
  cellmeta <- x@cellmeta
  if (!split.by %in% colnames(cellmeta)) {
    stop(paste(split.by, "is not found in 'cellmeta'"))
  }
  if (!is.factor(cellmeta[[split.by]]) && !is.character(cellmeta[[split.by]])) {
    stop(paste(split.by, "field is not a factor or character."))
  }
  factor.levels <- unique(cellmeta[[split.by]])
  object.list <- lapply(factor.levels, function(xx) {
    subset2(x, cells = rownames(subset(cellmeta, get(split.by) == xx)))
  })
  names(object.list) <- factor.levels
  return(object.list)
}

setMethod(
  f = "show",
  signature = "CellProject",
  definition = function(object) {
    cat("An object of class", class(x = object), "\n")
    cat("@data: ", ncol(object@data), "features across", nrow(object@data), "cells.\n")
    cat("@embeddings: ", "(", nrow(object@embeddings), ",", ncol(object@embeddings), ").\n")
    cat("@refined.embeddings: ", "(", nrow(object@refined.embeddings), ",", ncol(object@refined.embeddings), ").\n")
    cat("@cellmeta: ", "(", nrow(object@cellmeta), ",", ncol(object@cellmeta), ").\n")
    nn.dim <- dim(object@neighbors$nn.idx)
    if (is.null(nn.dim)) {
      cat("@neighbors: NULL\n")
    } else {
      cat("@neighbors: ", nn.dim[2], "nearest neighbors graph on", nn.dim[1], "cells.\n")
    }
  }
)

