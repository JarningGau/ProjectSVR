#' @include projection.R
NULL

#' A wrapper for mapping query data onto reference.
#'
#' @param seu.q A query Seurat object.
#' @param reference Reference model
#' @param assay.q The assay used for reference mapping.
#' @param ncores Number of threads for calculation.
#' @return A Seurat object. The projected reference embeddings were saved in a
#' dimension reduction object named 'ref.umap'. The gene set score were saved in
#' a new assay named 'UCell'.
#'
#' @export
#'
MapQuery <- function(seu.q, reference, assay.q = "RNA", ncores = 1) {
  ## check parameters
  if (!inherits(seu.q, "Seurat")) stop("seu.q argument must be a Seurat object")
  if (!is.list(reference) || !all(c("genes", "models") %in% names(reference))) {
    stop("reference argument must be a list containing 'genes' and 'models' elements")
  }
  if (!is.character(assay.q) || !(assay.q %in% names(seu.q@assays))) {
    stop("assay.q argument must be a character string specifying an existing assay in the 'seu.q' Seurat object")
  }
  if (!is.numeric(ncores) || ncores < 1) stop("ncores argument must be a positive integer")
  if (!all(c("bg.genes", "gene.sets") %in% names(reference$genes))) {
    stop("bg.genes or gene.sets not found in reference$genes element")
  }
  gene.sets <- reference$genes$gene.sets
  if (length(names(gene.sets)) == 0) stop("gene.sets names are empty")
  max_cores <- parallel::detectCores() / 2
  ncores <- min(ncores, max_cores)

  message("#### Compute Gene Set Score Matrix ####")
  bg.genes <- reference$genes$bg.genes
  genes.use <- intersect(rownames(seu.q), bg.genes)
  matched.ratio <- length(genes.use) / length(bg.genes)
  matched.ratio <- round(matched.ratio, 4)*100
  if (matched.ratio < 50) {
    stop(sprintf("Too less background genes (%s%%) were matched between query and reference. Please check the gene names in 'seu.q'", matched.ratio))
  }
  if (matched.ratio < 70) {
    warning(sprintf("Only %s%% background genes were matched between query and reference", matched.ratio))
  }
  if (matched.ratio) {
    message(sprintf("Only %s%% background genes were matched between query and reference", matched.ratio))
  }
  counts <- seu.q[[assay.q]]@counts[genes.use, ]
  if (!rlang::is_installed("UCell")) {
    stop("Please install UCell package (https://github.com/carmonalab/UCell).")
  }
  gss.mat <- UCell::ScoreSignatures_UCell(counts, features = gene.sets, ncores = ncores)
  colnames(gss.mat) <- paste0(names(gene.sets), "_UCell")

  message("#### Map Query to Reference ####")
  proj.obj <- ProjectNewdata(feature.mat = as.data.frame(gss.mat),
                             model = reference$models$umap,
                             do.norm = "L2", cores = ncores)
  seu.q[["UCell"]] <- Seurat::CreateAssayObject(data = gss.mat)
  seu.q[["ref.umap"]] <- Seurat::CreateDimReducObject(proj.obj@embeddings, key = "refUMAP_", assay = assay.q)
  return(seu.q)
}
