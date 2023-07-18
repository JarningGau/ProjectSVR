#' @include projection.R
#' @include label_transfer.R
NULL

#' A wrapper for mapping query data onto reference.
#'
#' @param seu.q A query Seurat object.
#' @param reference Reference model.
#' @param gss.method Method for calculation of gene set score. 'AUCell' or 'UCell'
#' are supported. Defalut: UCell
#' @param assay.q The assay used for reference mapping. Default: RNA
#' @param add.map.qual Whether add mapping quality metric (mean kNN distance). Default: FALSE
#' @param ncores Number of threads for calculation. Default: 1
#' @return A Seurat object. The projected reference embeddings were saved in a
#' dimension reduction object named 'ref.umap'. The gene set score were saved in
#' a new assay named 'SignatureScore'.
#' @concept reference_mapping
#' @export
#'
MapQuery <- function(seu.q, reference, gss.method = "UCell", assay.q = "RNA", add.map.qual = FALSE, ncores = 1) {
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
  seu.q <- ComputeModuleScore(seu.q, gene.sets, bg.genes, method = gss.method, cores = ncores)
  Seurat::DefaultAssay(seu.q) <- "SignatureScore"
  gss.mat <- Seurat::FetchData(seu.q, vars = rownames(seu.q))

  message("#### Map Query to Reference ####")
  proj.obj <- ProjectNewdata(feature.mat = gss.mat,
                             model = reference$models$umap,
                             do.norm = "L2", cores = ncores)
  seu.q[["ref.umap"]] <- Seurat::CreateDimReducObject(proj.obj@embeddings, key = "refUMAP_", assay = assay.q)
  if (add.map.qual) {
    message("#### Map Quality ####")
    proj.obj <- AddProjQual(object = proj.obj)
    seu.q$mean.knn.dist <- proj.obj@cellmeta$mean.knn.dist
    seu.q$mapQ.p.val <- proj.obj@cellmeta$p.val
    seu.q$mapQ.p.adj <- proj.obj@cellmeta$p.adj
  }
  return(seu.q)
}

#' A wrapper for transferring reference cell labels to the query via KNN method.
#'
#' @param seu.q A query Seurat object.
#' @param reference Reference model.
#' @param reduction.q Reduction name storing projected reference embeddings. Default is 'ref.umap'.
#' @param ref.emb.col The column names storing the reference embeddings in reference model.
#' @param ref.label.col The column name storing the cell labels to be transferred.
#' @param k The K param for KNN model. Default is 10.
#'
#' @return A Seurat object. The predicted cell type is store in 'knn.pred.celltype' column.
#' @concept label_transfer
#' @export
#'
LabelTransfer <- function(seu.q, reference, reduction.q = "ref.umap",
                          ref.emb.col = paste0("UMAP_", 1:2),
                          ref.label.col = "cluster.name", k = 10) {
  # 检查输入的 seu.q 是否为 Seurat 对象
  if(!inherits(seu.q, "Seurat")){
    stop("Input 'seu.q' must be a Seurat object.")
  }
  if (!(reduction.q %in% names(seu.q@reductions))) {
    stop("Input 'reduction.q' does not exist in seu.q.")
  }

  # 检查输入的 reference 是否包含指定的元素
  required_elements <- c("ref.cellmeta")
  missing_elements <- setdiff(required_elements, names(reference))
  if(length(missing_elements) > 0){
    stop(paste0("Input 'reference' is missing the following elements: ",
                paste(missing_elements, collapse = ", ")))
  }

  # 检查输入的 ref.emb.col 的长度是否合法
  if(length(ref.emb.col) != 2){
    stop("Input 'ref.emb.col' must contain two column names for embedding coordinates.")
  }

  # 检查输入的 ref.label.col 是否在 reference$ref.cellmeta$meta.data 中存在
  if(!(ref.label.col %in% colnames(reference$ref.cellmeta$meta.data))){
    stop("Input 'ref.label.col' does not exist in reference$ref.cellmeta.")
  }

  # 检查输入的 k 是否为正整数
  if(!is.numeric(k) || !k %% 1 == 0 || k <= 0){
    stop("Input 'k' must be a positive integer.")
  }
  query.emb <- seu.q[[reduction.q]]@cell.embeddings
  ref.cellmeta <- reference$ref.cellmeta$meta.data
  ref.emb <- ref.cellmeta[, ref.emb.col]
  ref.labels <- ref.cellmeta[[ref.label.col]]
  names(ref.labels) <- rownames(ref.cellmeta)
  knn.pred <- KnnLabelTransfer(query.emb = query.emb, ref.emb = ref.emb, ref.labels = ref.labels, k = k)
  seu.q$knn.pred.celltype <- knn.pred$labels
  seu.q$knn.pred.votes <- knn.pred$votes
  seu.q$knn.pred.perc <- knn.pred$perc
  ref.celltype.levels <- levels(ref.cellmeta[[ref.label.col]])
  if (!is.null(ref.celltype.levels)) {
    seu.q$knn.pred.celltype <- factor(seu.q$knn.pred.celltype, levels = ref.celltype.levels)
  }
  return(seu.q)
}
