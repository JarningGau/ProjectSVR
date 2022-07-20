#' @include utils.R
#' @include objects.R
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Multi Classifer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Balance sampling:
## https://github.com/Teichlab/celltypist/blob/main/celltypist/train.py

#' Training an ensemble classifier for cell type prediction.
#'
#' @import mlr3verse
#' @importFrom mlr3 as_task_classif lrn
#' @param feature.mat A signature score matrix, rows are cells, columns are features.
#' @param cell.types A data.frame recording cell types at different granularity.
#' @param do.norm Whether normalize the feature matrix. L1, L2, NULL. Default: NULL.
#' @param mlr3.model classif.svm, classif.xgboost, classif.randomForest.  Default: classif.svm
#' @param batch.size The number of cells for each model. Default: 5000
#' @param n.models The number of SVM model. Default: 100
#' @param balance.cell.type A boolen determines whether performing balance sampling. Default: True
#' @param cores The number of CPU for training. Default: -1, use all available threads.
#' @return a list of trained learners.
#' @export
FitEnsemblMultiClassif <- function(feature.mat, cell.types, do.norm = NULL, mlr3.model = "classif.svm", batch.size=5000,
                                   n.models=100, balance.cell.type = TRUE, cores=-1){
  ## check parameters
  if (!is.data.frame(feature.mat)) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  if (!all(rownames(feature.mat) == rownames(cell.types))) {
    stop("The rownames of 'feature.mat' and 'cell.types' are not matched.")
  }
  if (batch.size >= nrow(feature.mat)) {
    stop("Batch size is larger than the given cells. A smaller 'batch.size' is expected.")
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## check data sets: colnames should not contain "-"
  bad.cols <- grepl("-", colnames(feature.mat))
  if (any(bad.cols)) {
    warning("Bad colnames in 'feature.mat'. Change the '-' to '_'.")
    colnames(feature.mat) <- gsub("-", "_", colnames(feature.mat))
  }
  bad.cols <- grepl("-", colnames(cell.types))
  if (any(bad.cols)) {
    warning("Bad colnames in 'cell.types'. Change the '-' to '_'.")
    colnames(cell.types) <- gsub("-", "_", colnames(cell.types))
  }
  for (i in seq_along(colnames(cell.types))) {
    cell.types[[i]] %<>% as.factor()
  }
  ## 0. Normalize: by cells (rows)
  if (do.norm == "L1") {
    message("L1 normalization ...")
    feature.mat <- L1norm(feature.mat)
  } else if (do.norm == "L2") {
    message("L2 normalization ...")
    feature.mat <- L2norm(feature.mat)
  } else if (is.null(do.norm)) {
    message("Skip normalization ...")
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  ## 1. Create multi-classify tasks
  ## 1.1 sampling
  message(sprintf("Creating tasks: size=%s, n=%s", batch.size, n=n.models))
  .getP <- function(labels) {
    labels.freq <- as.vector(table(labels))
    names(labels.freq) <- names(table(labels))
    1 / (labels.freq[labels] * length(labels.freq))
  }
  ## 1.2 build tasks for each targets
  targets <- colnames(cell.types)
  tasks <- lapply(targets, function(ii) {
    if (balance.cell.type) {
      p <- .getP(cell.types[[ii]])
    } else {
      p <- NULL
    }
    pbapply::pblapply(1:n.models, function(jj) {
      row.ids <- sample(rownames(feature.mat), size = batch.size, replace = F, prob = p)
      data.use <- cbind(feature.mat[row.ids, , drop=FALSE], cell.types[row.ids, ii, drop=FALSE])
      as_task_classif(as.data.frame(data.use), target = ii)
    })
  })
  names(tasks) <- targets
  ## 2. Training
  message("Training ...")
  .train <- function(dataset, batch.id) {
    task = dataset[[batch.id]]
    if (!mlr3.model %in% mlr3::mlr_learners$keys()) {
      stop(sprintf("Undefined 'mlr3.model': %s.", mlr3.model))
    }
    learner = lrn(mlr3.model)
    learner$train(task)
    learner
  }
  celltype_lrns = lapply(tasks, function(celltype_i) {
    parallel::mclapply(1:n.models, function(j) .train(celltype_i, j), mc.cores=cores)
  })
  names(celltype_lrns) = names(tasks)
  model <- list(
    model = celltype_lrns,
    features = colnames(feature.mat)
  )
  class(model) <- "EmsemblMultiClassif"
  return(model)
}

#' Predict cell types
#'
#' @param feature.mat A data.frame containing signature scores, rows are cells, columns are features.
#' @param model The trained learners returned from `FitEnsemblMultiClassif()`.
#' @param do.norm Whether normalize the feature matrix. L1, L2, NULL. Default: NULL.
#' @param cores number of threads for prediction, -1 means all available threads. Default: -1
#' @return A data.frame contains predicted results.
#' @export
PredictNewdata <- function(feature.mat, model, do.norm=NULL, cores=-1){
  ## check parameters
  if (!is.data.frame(feature.mat)) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  bad.features <- setdiff(colnames(feature.mat), model$features)
  if (length(bad.features) > 0) {
    stop(sprintf("These features are not in the model: %s", paste(bad.features, collapse = " ")))
  }
  no.features <- setdiff(model$features, colnames(feature.mat))
  if (length(no.features) > 0) {
    warning(sprintf("These features are not in your inputs: %s\n. Filled with zeros.", paste(no.features, collapse = " ")))
    zeros <- matrix(rep(0, nrow(newdata)*length(no.features)), ncol = length(no.features)) %>% as.data.frame()
    colnames(zeros) <- no.features
    newdata <- cbind(feature.mat, zeros)[, model$features]
  } else {
    newdata <- feature.mat
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## Normalization: by cells (rows)
  if (do.norm == "L1") {
    message("L1 normalization ...")
    newdata <- L1norm(newdata)
  } else if (do.norm == "L2") {
    message("L2 normalization ...")
    newdata <- L2norm(newdata)
  } else if (is.null(do.norm)) {
    message("Skip normalization ...")
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  ## prediction
  n.models <- length(model$model[[1]])
  celltype_pred <- lapply(model$model, function(lrns) {
    celltype_i <- parallel::mclapply(
      X = 1:n.models,
      FUN = function(i) {
        lrns[[i]]$predict_newdata(newdata = as.data.frame(newdata))$data$response %>% as.character()
      }, mc.cores=cores)
    do.call(cbind, celltype_i)
  })
  names(celltype_pred) = names(model$model)
  ## integration
  int.fun <- function(xx) {
    celltype.counts <- sort(table(xx), decreasing = T)
    names(celltype.counts)[1]
  }
  celltype_inter <- lapply(names(celltype_pred), function(xx) {
    apply(celltype_pred[[xx]], 1, int.fun)
  })
  names(celltype_inter) <- names(celltype_pred)
  celltype_inter <- as.data.frame(celltype_inter)
  rownames(celltype_inter) <- rownames(feature.mat)
  return(celltype_inter)
}


# Over cluster via k-means for majority voting
OverCluster <- function(feature.mat, k){
  clusters <- stats::kmeans(feature.mat, centers = k)$cluster
  # clusters <- paste0("C", clusters)
  # clusters <- factor(clusters, levels = paste0("C", 1:k))
  return(clusters)
}


#' Majority voting to harmonise cell labels on given over clusters
#'
#' @param feature.mat A signature score matrix for over clustering, rows are cells, columns are features.
#' @param over.clusters A named vector contains clusters for majority vote. Default: NULL
#' @param cell.types A data.frame recording cell types at different granularity.
#' @param k k of kmeans for over cluster. Default: 20
#' @param min.prop The minimum proportion of cells required to support naming of
#' the subcluster by a cell type. Default: 0
#' @return A data.frame contains consensus predicted cell types.
#' @export
MajorityVote <- function(feature.mat = NULL, over.clusters = NULL, cell.types, k = 20, min.prop = 0){
  ## check parameters
  if (is.null(feature.mat) && is.null(over.clusters)) {
    stop("Must provide one of the 'feature.mat' and 'over.clusters'.")
  }
  if (!is.null(over.clusters) && is.null(names(over.clusters))) {
    stop("'over.clusters' should be named.")
  }
  if (!is.null(feature.mat) && is.null(rownames(feature.mat))) {
    stop("No rownames of 'feature.mat'")
  }
  ## Over cluster
  if (is.null(over.clusters)) {
    common.cells <- intersect(rownames(feature.mat), rownames(cell.types))
    if (length(common.cells) != nrow(cell.types)) {
      stop("Inconsistent rownames between 'feature.mat' and 'cell.types'.")
    }
    over.clusters <- OverCluster(feature.mat, k = k)
  } else {
    common.cells <- intersect(names(over.clusters), rownames(cell.types))
    if (length(common.cells) != nrow(cell.types)) {
      stop("Inconsistent cell names between 'over.clusters' and 'cell.types'.")
    }
  }
  ## Majority vote
  all.vars <- colnames(cell.types)
  cluster.levels <- unique(over.clusters)
  major.votes <- lapply(cluster.levels, function(clu) {
    cells.use <- names(over.clusters[over.clusters == clu])
    row.mv <- sapply(seq_along(all.vars), function(ii) {
      freq <- table(cell.types[cells.use, all.vars[ii]]) / length(cells.use)
      top1 <- sort(freq, decreasing = TRUE)[1] %>% names()
      ifelse(top1 > min.prop, top1, "Heterogeneous") ## -> return
    })
    names(row.mv) <- paste0(all.vars, ".major_votes")
    row.mv
  })
  major.votes <- do.call(rbind, major.votes)
  rownames(major.votes) <- cluster.levels # rows = clusters, columns = xx.major_votes
  ## Mapping clusters to major cell types
  cell.types$over.clusters <- as.character(over.clusters[rownames(cell.types)])
  pred.df <- major.votes[cell.types$over.clusters, , drop=FALSE]
  rownames(pred.df) <- rownames(cell.types)
  return(as.data.frame(pred.df))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### KNN based label transfer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Label transfer via KNN
#'
#' @param pred.emb The predicted embedding matrix of query data.
#' @param ref.emb The embedding matrix of reference atlas.
#' @param ref.celltype A named vector <name, value> = <cellname, celltype>.
#' @param k K nearest neighbors used. Default: 100
#'
#' @return A data.frame contains prediction results.
#'
#' @importFrom magrittr %>%
#' @export
#'
KnnLabelTransfer <- function(pred.emb, ref.emb, ref.celltype, k=100) {
  ## check parameters
  if (is.null(names(ref.celltype))) {
    stop("'ref.celltype' should be named.")
  }
  common.cells <- intersect(rownames(ref.emb), names(ref.celltype))
  if (length(common.cells) != nrow(ref.emb) || length(common.cells) != length(ref.celltype)) {
    stop("The rownames of 'ref.emb' should be consistant with the names of 'ref.celltype'.")
  }
  # reorder ref.celltype
  ref.celltype <- ref.celltype[rownames(ref.emb)]
  ## 1. build kNN graph on ref.emb
  ## 2. find k nearest neightbors for query on ref.emb
  nn.ranked <- NNHelper(data = ref.emb, query = pred.emb, k = k+1, method = "rann")
  ## 3. label transfer
  results <- pbapply::pblapply(1:nrow(nn.ranked$nn.idx), function(i) {
    xx <- nn.ranked$nn.idx[i, ]
    n.celltype <- sort(table(ref.celltype[xx]), decreasing = TRUE)
    most.voted <- n.celltype[1]
    data.frame(
      celltype.pred = names(most.voted),
      votes = most.voted,
      perc = most.voted / length(xx)
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- nn.ranked$query.cell.names
  return(results)
}
