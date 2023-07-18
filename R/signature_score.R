#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Compute Module Score ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate gene module scores
#'
#' @param x gene expression matrix, rows are genes, columns are cells. Can be any
#' format, UMI, CPM, TPM, etc.
#' @param ... Arguments passed to other methods.
#' @return A signature score matrix or Seurat object.
#' @concept compute_module_score
#' @export
ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')


#' @param gene.sets a list of gene sets in data.frame or named list.
#' @param bg.genes background genes for calculating signature score, NULL means use
#' all genes. Default: NULL
#' @param method method for calculating signature score, support "AUCell" or "UCell".
#' Default: UCell
#' @param min.size The minimal genes of the gene sets. The size of gene sets less
#' than this value were ignored.  Default: 20
#' @param batch.size The number of cells were calculated for each batch.
#' This parameter is for the parallel calculation in adopted in 'AUCell' method.
#' Default: 500
#' @param cores number of threads for parallel computing. Default: 1
#'
#' @rdname ComputeModuleScore
#' @concept compute_module_score
#' @export
ComputeModuleScore.default <- function(x, gene.sets, bg.genes=NULL, method="UCell",
                                       min.size=20, batch.size=500, cores=1, ...) {
  # Check if UCell and AUCell packages are installed
  if (method == "AUCell" && !requireNamespace("AUCell", quietly = TRUE)) {
    stop("AUCell package is not installed. Please install the package using BiocManager::install('AUCell').")
  } else if (method == "UCell" && !requireNamespace("UCell", quietly = TRUE)) {
    stop("UCell package is not installed. Please install the package using remotes::install_github('carmonalab/UCell', ref='v1.3').")
  }
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  # filter the gene sets
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  # check the bg.genes
  if (!is.null(bg.genes)) {
    genes.use <- intersect(bg.genes, rownames(x))
    matched.ratio <- length(genes.use) / length(bg.genes)
    matched.ratio <- round(matched.ratio, 4)*100
    if (matched.ratio < 50) {
      stop(sprintf("Too less background genes (%s%%) were matched between query and reference. Please check the gene names in 'x'", matched.ratio))
    }
    if (matched.ratio < 70) {
      warning(sprintf("Only %s%% background genes were matched between query and reference", matched.ratio))
    } else {
      message(sprintf("%s%% background genes were matched between query and reference", matched.ratio))
    }
    x <- x[genes.use, ]
  }
  # calculate signature scores
  if (method == "AUCell") {
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
    auc_scores <- do.call(cbind, auc_scores)
    auc_scores <- as.data.frame(t(auc_scores))
  } else if (method == "UCell") {
    auc_scores <- UCell::ScoreSignatures_UCell(x, features = gene.sets, ncores = cores)
  } else {
    stop(sprintf("%s is not supported! Please set method to 'AUCell' or 'UCell'.", method))
  }
  colnames(auc_scores) <- names(gene.sets)
  return(auc_scores)
}


#' @rdname ComputeModuleScore
#' @param assay Name of the seurat object assay.
#' @concept compute_module_score
#' @export
ComputeModuleScore.Seurat <- function(x, gene.sets, bg.genes=NULL, method="UCell",
                                      min.size=20, batch.size=500, cores=1,
                                      assay = Seurat::DefaultAssay(x), ...) {
  dge <- x[[assay]]@counts
  auc_scores <- ComputeModuleScore.default(x = dge, gene.sets, bg.genes, method,
                                           min.size, batch.size, cores)
  x[["SignatureScore"]] <- Seurat::CreateAssayObject(data = t(auc_scores))
  return(x)
}


