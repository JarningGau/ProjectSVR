# Consensus NMF analysis via python package `cNMF` through `reticulate`.
#
#' @importFrom utils head write.table
#' @importFrom magrittr "%>%"
NULL

#' Find optimal components or topics for cNMF
#'
#' @param counts.fn path to the cell x gene counts file. This is expected to be
#' a tab-delimited text file or a \code{anndata} object saved in the h5ad format.
#' @param components vector, list of K values that will be tested for cNMF.
#' @param out.path the output directory into which all results will be placed.
#' Default: current path.
#' @param run.name a subdirectory out.path/run.name will be created and all output
#' files will have name as their prefix. Default: current timestamp.
#' @param n.iter number of NMF iterations to run for each K. Default: 100.
#' @param tpm.fn If provided, load tpm data from file. Otherwise will compute it
#' from the counts file. Default: NULL
#' @param n.var.genes (optional) the number of highest variance genes that will
#' be used for running the factorization. Default: 2000
#' @param genes.fn (optional)  List of over-dispersed genes to be used for the
#' factorization steps. One gene per line. If not provided, over-dispersed
#' genes will be calculated automatically and the number of genes to use can be
#' set by the `n.var.genes` parameter below. Default: NULL
#' @param seed the master seed that will be used to generate the individual seed
#' for each NMF replicate. Default: 1024
#' @param cores specifies how many cores can be used in parallel. Default: all
#' available cores detected by `parallel::detectCores()`.
#' @export
FindOptimalK <- function(counts.fn, components, tpm.fn=NULL, out.path=NULL, run.name=NULL,
                         n.iter=100, n.var.genes=2000, genes.fn=NULL, seed=1024, cores=-1) {
  # check parameters
  if (!file.exists(counts.fn)) {
    stop(sprintf("File %s does not exist.", counts.fn))
  }
  if (is.null(components)) {
    stop("place provide @param components.")
  }
  if (!is.null(genes.fn) && !file.exists(genes.fn)) {
    stop(sprintf("File %s does not exist.", genes.fn))
  }
  if (!is.null(tpm.fn) && !file.exists(tpm.fn)) {
    stop(sprintf("File %s does not exist.", tpm.fn))
  }
  out.path <- ifelse(is.null(out.path), file.path(getwd(), "cnmf_tmp"), out.path)
  run.name <- ifelse(is.null(run.name), as.numeric(Sys.time()), run.name)
  run.name <- as.character(run.name)
  if (!dir.exists(out.path)) {
    dir.create(out.path, recursive = TRUE)
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())

  cnmf <- reticulate::import("cnmf")
  # init cNMF object
  cnmf_obj = cnmf$cNMF(output_dir=out.path, name=run.name)
  # prepare step:
  # 1. select high variable genes if genes.fn is null.
  # 2. normalization via variance scale.

  components %<>% as.integer()
  n.iter %<>% as.integer()
  seed %<>% as.integer()
  if (is.null(genes.fn)) {
    n.var.genes %<>% as.integer()
    cnmf_obj$prepare(counts_fn=counts.fn, components=components, n_iter=n.iter, tpm_fn=tpm.fn, seed=seed, num_highvar_genes = n.var.genes)
  } else {
    cnmf_obj$prepare(counts_fn=counts.fn, components=components, n_iter=n.iter, tpm_fn=tpm.fn, seed=seed, genes_file = genes.fn)
  }

  # factorization
  cores %<>% as.integer()
  parallel::mclapply(
    X = 0:(cores-1),
    FUN = function(i) {
      i %<>% as.integer()
      cnmf_obj$factorize(worker_i=i, total_workers=cores)
      invisible()
    },
    mc.cores = cores
  )
  # merge replicates
  cnmf_obj$combine()
  # plot
  cnmf_obj$k_selection_plot()

  invisible()
}


#' Run consensus NMF
#'
#' @param counts.fn path to the cell x gene counts file. This is expected to be
#' a tab-delimited text file or a `anndata` object saved in the h5ad format.
#' @param K int, number of components (topics or dimensions) for cNMF.
#' @param out.path the output directory into which all results will be placed.
#' Default: current path.
#' @param run.name a subdirectory out.path/run.name will be created and all output
#' files will have name as their prefix. Default: current timestamp.
#' @param n.iter number of NMF iterations to run for each K. Default: 100.
#' @param n.var.genes (optional) the number of highest variance genes that will
#' be used for running the factorization. Default: 2000
#' @param genes.fn (optional)  List of over-dispersed genes to be used for the
#' factorization steps. One gene per line. If not provided, over-dispersed
#' genes will be calculated automatically and the number of genes to use can be
#' set by the `n.var.genes` parameter below. Default: NULL
#' @param seed the master seed that will be used to generate the individual seed
#' for each NMF replicate. Default: 1024
#' @param cores specifies how many cores can be used in parallel. Default: all
#' available cores detected by `parallel::detectCores()`.
#' @param n.top.genes number of the genes with the highest loadings in each gene
#' expression program (GEP). Default: 100
#' @param local.density.cutoff the threshold on average distance to K nearest
#' neighbors to use. 2.0 or above means that nothing will be filtered out. Default: 0.5
#' @param local.neighborhood.size Percentage of replicates to consider as nearest
#' neighbors for local density filtering. E.g. if you run 100 replicates, and set
#' this to 0.3, 30 nearest neighbors will be used for outlier detection. Default: 0.3
#' @param show.clustering whether or not the clustergram image is output. Default: FALSE
#' @export
RunCNMF <- function(counts.fn, K, out.path=NULL, run.name=NULL, n.iter=100,
                    n.var.genes=2000, genes.fn=NULL, seed=1024, cores=-1, n.top.genes=100,
                    local.density.cutoff=0.5, local.neighborhood.size=0.3,
                    show.clustering=FALSE) {
  # check parameters
  if (!file.exists(counts.fn)) {
    stop(sprintf("File %s does not exist.", counts.fn))
  }
  if (is.null(K)) {
    stop("place provide @param components.")
  }
  if (!is.null(genes.fn) && !file.exists(genes.fn)) {
    stop(sprintf("File %s does not exist.", genes.fn))
  }
  out.path <- ifelse(is.null(out.path), file.path(getwd(), "cnmf_tmp"), out.path)
  run.name <- ifelse(is.null(run.name), as.numeric(Sys.time()), run.name)
  run.name <- as.character(run.name)
  if (!dir.exists(out.path)) {
    dir.create(out.path, recursive = TRUE)
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())

  cnmf <- reticulate::import("cnmf")
  # init cNMF object
  cnmf_obj = cnmf$cNMF(output_dir=out.path, name=run.name)
  # prepare:
  # 1. select high variable genes if genes.fn is null.
  # 2. normalization via variance scale.
  K %<>% as.integer()
  n.iter %<>% as.integer()
  seed %<>% as.integer()
  if (is.null(genes.fn)) {
    n.var.genes %<>% as.integer()
    cnmf_obj$prepare(counts_fn=counts.fn, components=K, n_iter=n.iter, seed=seed, num_highvar_genes = n.var.genes)
  } else {
    cnmf_obj$prepare(counts_fn=counts.fn, components=K, n_iter=n.iter, seed=seed, genes_file = genes.fn)
  }
  # factorization
  parallel::mclapply(
    X = 0:(cores-1),
    FUN = function(i) {
      i %<>% as.integer()
      cnmf_obj$factorize(worker_i=i, total_workers=cores)
      invisible()
    },
    mc.cores = cores
  )
  # merge replicates
  cnmf_obj$combine()
  # consensus NMF
  cnmf_obj$consensus(k=K, density_threshold=local.density.cutoff,
                     local_neighborhood_size=local.neighborhood.size,
                     show_clustering=show.clustering)
  # extract genes with high loadings in each GEP
  local.density.cutoff <- sub("\\.", "_", local.density.cutoff)
  gene.loading.fn <- sprintf("%s.gene_spectra_score.k_%s.dt_%s.txt", run.name, K, local.density.cutoff)
  gene.loading.df <- data.table::fread(file.path(out.path, run.name, gene.loading.fn), header = TRUE, sep = "\t")
  gene.loading.df <- gene.loading.df[,-1]
  top.genes <- apply(gene.loading.df, 1, function(yy) head(names(sort(yy, decreasing = T)), n.top.genes) )
  colnames(top.genes) <- paste0("program_", 1:K)
  top.genes.fn <- sprintf("top%s_genes.k_%s.dt_%s.txt", n.top.genes, K, local.density.cutoff)
  write.table(top.genes, file = file.path(out.path, run.name, top.genes.fn), sep = "\t", quote = F, row.names = FALSE)

  invisible()
}
