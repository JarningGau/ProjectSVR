#' @include utils.R
#' @include objects.R
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Regression Model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Training an ensemble SVM model for embedding coordinates regression.
#'
#' @import mlr3verse
#' @importFrom mlr3 as_task_regr lrn
#'
#' @param feature.mat A signature score matrix, rows are cells, columns are features.
#' @param emb.mat The embedding matrix, rows are cells, columns are dimensions.
#' @param cell.types A named vector recording cell types.
#' @param do.norm Whether normalize the feature matrix. L1, L2, NULL. Default: 'L2'
#' @param batch.size The number of cells for each model. Default: 5000
#' @param n.models The number of SVM model. Default: 100
#' @param balance.cell.type A boolen determines whether performing balance sampling. Default: FALSE
#' @param cores The number of CPU for training. Default: -1, use all available threads.
#' @return a list of trained learners.
#' @concept training_model
#' @export
FitEnsembleSVM <- function(feature.mat, emb.mat, cell.types=NULL, do.norm='L2', batch.size=5000, n.models=100, balance.cell.type=FALSE, cores=-1) {
  ## check parameters
  if (!is.data.frame(feature.mat)) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  outliers.1 <- setdiff(rownames(feature.mat), rownames(emb.mat))
  outliers.2 <- setdiff(rownames(emb.mat), rownames(feature.mat))
  if (length(outliers.1)) {
    stop("Some cells in feature.mat are not in emb.mat.")
  }
  if (length(outliers.2)) {
    stop("Some cells in emb.mat are not in feature.mat.")
  }
  if (batch.size >= nrow(feature.mat)) {
    stop("Batch size is larger than the samples. A smaller 'batch.size' is expected.")
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## check data sets: colnames should not contain "-"
  bad.cols <- grepl("-", colnames(feature.mat))
  if (any(bad.cols)) {
    warning("Bad colnames in 'feature.mat'. Change the '-' to '_'.")
    colnames(feature.mat) <- gsub("-", "_", colnames(feature.mat))
  }
  bad.cols <- grepl("-", colnames(emb.mat))
  if (any(bad.cols)) {
    warning("Bad colnames in 'emb.mat'. Change the '-' to '_'.")
    colnames(emb.mat) <- gsub("-", "_", colnames(emb.mat))
  }
  ## 0. Normalize: by cells (rows)
  if (is.null(do.norm)) {
    message("Skip normalization ...")
  } else if (do.norm == "L1") {
    message("L1 normalization ...")
    feature.mat <- L1norm(feature.mat)
  } else if (do.norm == "L2") {
    message("L2 normalization ...")
    feature.mat <- L2norm(feature.mat)
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  ## 1. Create regression tasks
  .getP <- function(labels) {
    labels.freq <- as.vector(table(labels))
    names(labels.freq) <- names(table(labels))
    1 / (labels.freq[labels] * length(labels.freq))
  }
  ## 1.1 sampling
  message(sprintf("Creating regression tasks: size=%s, n=%s", batch.size, n=n.models))
  if (balance.cell.type & !is.null(cell.types)) {
    message("Balanced sampling.")
    p <- .getP(cell.types)
  } else {
    p <- NULL
  }
  train.datasets <- pbapply::pblapply(1:n.models, function(i) {
    row.ids <- sample(rownames(feature.mat), size = batch.size, replace = F, prob = p)
    as.data.frame(cbind(feature.mat[row.ids, , drop=FALSE], emb.mat[row.ids, , drop=FALSE]))
  })
  ## 1.2 build tasks for each targets
  vars <- colnames(feature.mat)
  targets <- colnames(emb.mat)
  tasks <- lapply(targets, function(coord_i) {
    pbapply::pblapply(train.datasets, function(data_j) {
      as_task_regr(data_j[, c(vars, coord_i)], target = coord_i)
    })
  })
  names(tasks) <- targets
  ## 2. Training
  .train <- function(dataset, batch.id) {
    task = dataset[[batch.id]]
    learner = lrn("regr.svm")
    learner$train(task)
    learner
  }
  message("Regression ...")
  coords_lrns = lapply(tasks, function(coord_i) {
    parallel::mclapply(1:n.models, function(j) .train(coord_i, j), mc.cores=cores)
  })
  names(coords_lrns) = names(tasks)
  model <- list(
    model = coords_lrns,
    features = colnames(feature.mat)
  )
  class(model) <- "Regression"
  return(model)
}


#' Predicting the embedding coordinates via the trained SVM model.
#'
#' @param feature.mat A data.frame containing signature scores, rows are cells, columns are features.
#' @param model The trained learners returned from `FitEnsembleSVM()`.
#' @param do.norm Whether normalize the feature matrix. L1, L2, NULL. Default: 'L2'
#' @param int.fun The function for integration of predicted embedding by different learners:
#' mean, median, or any self defined function given a vector and returns a value. Default: median.
#' @param cores number of threads for prediction, -1 means all available threads. Default: -1
#' @return \code{CellProject} object containing source data and predicted embedding matrix.
#' @concept reference_mapping
#' @export
ProjectNewdata <- function(feature.mat, model, do.norm='L2', int.fun=stats::median, cores=-1) {
  ## check parameters
  if (!is.data.frame(feature.mat)) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  bad.features <- setdiff(colnames(feature.mat), model$features)
  if (length(bad.features) > 0) {
    warning(sprintf("These features are not in the model: %s. Dropping them (it).", paste(bad.features, collapse = " ")))
    feature.mat <- feature.mat[, !(colnames(feature.mat) %in% bad.features), drop = FALSE]
  }
  no.features <- setdiff(model$features, colnames(feature.mat))
  if (length(no.features) > 0) {
    warning(sprintf("These features are not in your inputs: %s\n. Filling with zeros.", paste(no.features, collapse = " ")))
    zeros <- matrix(rep(0, nrow(newdata)*length(no.features)), ncol = length(no.features)) %>% as.data.frame()
    colnames(zeros) <- no.features
    newdata <- cbind(feature.mat, zeros)[, model$features]
  } else {
    newdata <- feature.mat
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## Normalization: by cells (rows)
  if (is.null(do.norm)) {
    message("Skip normalization ...")
  } else if (do.norm == "L1") {
    message("L1 normalization ...")
    newdata <- L1norm(newdata)
  } else if (do.norm == "L2") {
    message("L2 normalization ...")
    newdata <- L2norm(newdata)
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  ## prediction
  n.models <- length(model$model[[1]])
  coords_pred <- lapply(model$model, function(lrns) {
    coords_i <- parallel::mclapply(
      X = 1:n.models,
      FUN = function(i) {
        lrns[[i]]$predict_newdata(newdata = as.data.frame(newdata))$data$response
      }, mc.cores=cores)
    do.call(cbind, coords_i)
  })
  names(coords_pred) = names(model$model)
  ## integration
  coords_inter <- lapply(names(coords_pred), function(xx) {
    apply(coords_pred[[xx]], 1, int.fun)
  })
  names(coords_inter) <- names(coords_pred)
  coords_inter <- as.data.frame(coords_inter)
  rownames(coords_inter) <- rownames(feature.mat)
  ## returns
  newCellProject(embeddings = coords_inter, data = feature.mat)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Mapping Quanlity ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add projection quality scores for the query cells.
#'
#' @param object A \code{CellProject} object
#' @param k K nearest neighbors used. Default: 20
#' @param repeats The number of sample times for generate background projection qualities.
#' @return A \code{CellProject} object
#' @concept reference_mapping
#' @references \url{https://www.nature.com/articles/s41467-021-25957-x}
#' @export
AddProjQual <- function(object, k=20, repeats=1e4) {
  pred.emb <- object@embeddings
  golden.emb <- object@data

  if (!all(rownames(golden.emb) == rownames(pred.emb))) {
    pred.emb <- pred.emb[rownames(golden.emb), ]
  }

  message("Building k-NN graph in feature space ...")
  nn.golden = NNHelper(golden.emb, k = k, method = "sklearn", metric = "cosine")

  message("Calculating euclidean metric on pred.emb ...")
  pred.dist <- as.matrix(stats::dist(pred.emb))

  message("Calculating cell-based mapping quanlity metrics ...")
  mean.knn.dist <- pbapply::pbsapply(1:nrow(pred.emb), function(i) {
    mean(pred.dist[i, nn.golden$nn.idx[i, ]])
  })
  names(mean.knn.dist) <- rownames(pred.emb)

  mean.krand.dist <- pbapply::pbsapply(1:repeats, function(i) {
    rand.cells <- sample(1:ncol(pred.dist), size = k+1, replace = F)
    mean(pred.dist[rand.cells[1], rand.cells[-1]])
  })
  p.val <- pbapply::pbsapply(mean.knn.dist, function(q) sum(q > mean.krand.dist) / length(mean.krand.dist) )
  p.adj <- stats::p.adjust(p.val, method = "BH")
  results <- list(
    mean.knn.dist = mean.knn.dist,
    # mean.krand.dist = mean.krand.dist, # this field is not equal to others.
    p.val = p.val,
    p.adj = p.adj
  )
  results <- as.data.frame(results)
  object@neighbors <- nn.golden
  cellmeta <- object@cellmeta
  ## check duplicated colnames in cellmeta
  res.names <- colnames(results)
  kept.names <- setdiff(colnames(cellmeta), res.names)
  cellmeta <- cellmeta[, kept.names, drop=FALSE]
  object@cellmeta <- cbind(cellmeta, results[rownames(cellmeta), ])
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Refine Mapping Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Refine low confidence mapping coords.
#'
#' @param object A \code{CellProject} object
#' @param p.val.cutoff Interpolation for the coordinates of cells with p.value
#' above this cutoff. Default: 1
#' @param p.adj.cutoff Interpolation for the coordinates of cells with p.value
#' above this cutoff. Default: 0.05
#' @param k K nearest neighbors for interpolation. Default: 10
#' @return A \code{CellProject} object
#' @concept reference_mapping
#' @export
RefineProjection <- function(object, p.val.cutoff = 1, p.adj.cutoff = 0.05, k = 10){
  if (is.null(object@neighbors)) {
    stop("No neighbors object provided.")
  }
  data <- object@data
  nn.dists <- object@neighbors$nn.dists
  query.cell.names <- object@neighbors$query.cell.names
  rownames(nn.dists) <- query.cell.names
  embeddings <- object@embeddings
  cellmeta <- object@cellmeta
  refined.embeddings <- embeddings
  if (k > ncol(nn.dists)) {
    warning("K =", k, "was too big, reset to", ncol(nn.dists), "according to the built knn graph.")
    k <- ncol(nn.dists)
  }
  sel.rows <- cellmeta$p.val > p.val.cutoff | cellmeta$p.adj > p.adj.cutoff
  low.qual.map <- cellmeta[sel.rows, ]
  # low.qual.map <- subset(cellmeta, p.val > p.val.cutoff | p.adj > p.adj.cutoff)
  if (nrow(low.qual.map)) {
    cat("There are", nrow(low.qual.map), "cell(s) with low quanlity mapping coordinates under p <=", p.val.cutoff,
        "and p.adj <=", p.adj.cutoff, "\n")
    ## 1) distance matrix
    low.qual.cells <- rownames(low.qual.map)
    other.cells <- setdiff(rownames(cellmeta), low.qual.cells)
    D <- cosine.dist(as.matrix(data))[low.qual.cells, other.cells]
    ## 2) sigma
    sigma <- nn.dists[low.qual.cells, k]
    sigma <- matrix(rep(sigma, each=length(other.cells)), ncol = length(other.cells), byrow = T)
    ## 3) affinity matrix
    A <- exp(-(D / sigma)^2)
    ## 4) normalize, rows are low.qual.cells, columns are other cells
    A <- t(apply(A, 1, function(xx) xx/sum(xx) ))
    ## 5) interpolation
    U <- A %*% embeddings[other.cells, ]
    refined.embeddings[rownames(U), ] <- U
  }
  object@refined.embeddings <- refined.embeddings
  return(object)
}

