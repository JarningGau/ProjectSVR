#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Compute Module Score ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate gene module scores
#'
#' @param x gene expression matrix, rows are genes, columns are cells. Can be any format, UMI, CPM, TPM, etc.
#' @param ... Arguments passed to other methods.
#' @return A signature score matrix or Seurat object.
#' @concept compute_module_score
#' @export
ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')


#' @param gene.sets a list of gene sets in data.frame or named list.
#' @param min.size The minimal genes of the gene sets. The size of gene sets less than this value were ignored.  Default: 20
#' @param batch.size The number of cells were calculated for each batch. Default: 500
#' @param cores number of threads for parallel computing. Default: 1
#'
#' @rdname ComputeModuleScore
#' @concept compute_module_score
#' @export
ComputeModuleScore.default <- function(x, gene.sets, min.size=20, batch.size=500, cores=1, ...) {
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  n.cells <- ncol(x)
  batches <- floor((1:n.cells-1) / batch.size)
  batch.levels <- unique(batches)
  aucell <- function(i) {
    dge.tmp <- x[, batches == i]
    cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores=1, plotStats=F, verbose = F)
    auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores=1, verbose = F)
    AUCell::getAUC(auc)
  }
  auc_scores <- parallel::mclapply(batch.levels, aucell, mc.cores = cores)
  do.call(cbind, auc_scores)
}


#' @rdname ComputeModuleScore
#' @param assay Name of the seurat object assay.
#' @concept compute_module_score
#' @export
ComputeModuleScore.Seurat <- function(x, gene.sets, min.size=20, batch.size=500, cores=1, assay = Seurat::DefaultAssay(x), ...) {
  dge <- x[[assay]]@counts
  ras_mat <- ComputeModuleScore.default(x = dge, gene.sets, min.size, batch.size, cores)
  x[["AUCell"]] <- Seurat::CreateAssayObject(data = ras_mat)
  return(x)
}


