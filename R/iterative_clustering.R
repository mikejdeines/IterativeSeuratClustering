initial_clustering <- function(seurat_object){
  #' Performs an initial clustering iteration on a Seurat object.
  #' Runs Leiden clustering at low resolution to find the initial clusters and runs one clustering iteration.
  #' @param seurat_object a normalized, integrated Seurat object
  #' @returns a Seurat object with initial clusters in the "leiden_clusters" slot
  require(Seurat)
  require(scCustomize)
  seurat_object <- FindClusters(seurat_object, resolution = 0.01, algorithm = 4, method = "igraph", verbose = FALSE)
  objs <- SplitObject(seurat_object, "seurat_clusters")
  objs <- objs[order(names(objs))]
  samples <- lapply(objs, leiden_clustering)
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}
clustering_iteration <- function(seurat_object, min_score, cluster_size, pct.1){
  #' Performs an iteration of Leiden clustering to a Seurat object.
  #' @param seurat_object a normalized, integrated Seurat object
  #' @param min_score minimum clustering score
  #' @param cluster_size minimum cluster size
  #' @param pct.1 fraction of gene expression in the overexpressing cluster
  #' @returns a Seurat object with clusters in the "leiden_clusters" slot
  require(Seurat)
  require(scCustomize)
  Idents(seurat_object) <- seurat_object$leiden_clusters
  objs <- SplitObject(seurat_object, "leiden_clusters")
  objs <- objs[order(names(objs))]
  samples <- lapply(objs, leiden_clustering, score_limit = min_score, min_size = cluster_size, pct.1 = pct.1)
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}
iterative_clustering <- function(seurat_object, max_iterations = 10, min_score = 150, cluster_size = 20, pct.1 = 0.5, dims_use = 1:30, reduction.name = "pca"){
  #' Performs iterative clustering on a Seurat object to find clusters expressing significant DE genes.
  #' @param seurat_object a normalized, integrated Seurat object
  #' @param max_iterations number of clustering iterations to perform
  #' @param min_score minimum cluster score
  #' @param cluster_size minimum cluster size
  #' @param pct.1 fraction of gene expression in the overexpressing cluster
  #' @param dims_use number of dimensions to use for neighbors graph
  #' @reduction.name name of the dimensional reduction used to find neighbors
  #' @returns a Seurat object with clusters in the "leiden_clusters" slot
  seurat_object <- initial_clustering(seurat_object)
  cluster_sizes <- data.frame()
  for (i in 2:max_iterations-1){
    old_clusters <- length(levels(seurat_object$leiden_clusters))
    seurat_object <- clustering_iteration(seurat_object, min_score, cluster_size, pct.1)
    seurat_object$leiden_clusters <- as.factor(seurat_object$leiden_clusters)
    num_clusters <- length(levels(seurat_object$leiden_clusters))
    if (old_clusters - num_clusters == 0) break()
  }
  return(seurat_object)
}
leiden_clustering <- function(seurat_object, num_clusters = 2, score_limit = 150, min_size = 20, pct.1 = 0.5, dims.use = dims_use, reduction = reduction.name){
  #' Leiden clustering into a set number of clusters. Varies resolution until the number of clusters is achieved.
  #' Requires a conda environment with igraph and leidenalg installed.
  #' @param seurat_object a normalized, integrated Seurat object
  #' @param num_clusters number of desired clusters. Default = 2.
  #' @param score_limit minimum score for a cluster. Default = 150.
  #' @param min_size minimum size of a cluster. Default = 20.
  #' @param pct.1 fraction of gene expression in the overexpressing cluster. Default = 0.5.
  #' @returns a Seurat object with clusters in the "leiden_clusters" slot
  require(Seurat)
  seurat_object$starting_clusters <- Idents(seurat_object)
  if (ncol(seurat_object) < 3) {
  message("Skipping: too few cells for clustering.")
  seurat_object$leiden_clusters <- seurat_object$starting_clusters
  return(seurat_object)
  }
  available_pcs <- ncol(Embeddings(seurat_object, reduction = reduction))
  dims_to_use <- dims.use[dims.use <= available_pcs]
  k_val <- min(20, cell_count - 1)
  seurat_object <- FindNeighbors(seurat_object, dims = dims_to_use, reduction = reduction, verbose = FALSE, k.param = k_val)
  initial_resolution = 1
  seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph", verbose = FALSE)
  cluster_count <- length(levels(Idents(seurat_object)))
  while (cluster_count != num_clusters){
    if (cluster_count > num_clusters){
      initial_resolution = initial_resolution/10
      seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph", verbose = FALSE)
      cluster_count <- length(levels(seurat_object$seurat_clusters))
    }
    else{
      initial_resolution = initial_resolution*5
      seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph", verbose = FALSE)
      cluster_count <- length(levels(seurat_object$seurat_clusters))
    }
  }
  seurat_object$leiden_clusters <- Idents(seurat_object)
  seurat_object$leiden_clusters <- paste(seurat_object$starting_clusters, "_", seurat_object$leiden_clusters, sep = "")
  score <- clustering_score(seurat_object, pct.1 = pct.1)
  if (score < score_limit){
    seurat_object$leiden_clusters <- seurat_object$starting_clusters
  }
  sizes <- data.frame(table(seurat_object$leiden_clusters))
  for (i in 1:length(rownames(sizes))){
    if (sizes$Freq[[i]] < min_size){
      seurat_object$leiden_clusters <- seurat_object$starting_clusters
    }
  }
  return(seurat_object)
}
clustering_score <- function(seurat_object, pct.1){
  #' Calculates cluster score based on DE genes. Determines DE genes between two clusters with a Wilcoxon signed rank test.
  #' Cluster score is the sum of -log10 of Bonferroni-corrected p values. The maximum score any gene can contribute is 20.
  #' @param seurat_object a normalized, integrated Seurat object
  #' @param pct.1 fraction of gene expression in the overexpressing cluster
  #' @returns a clustering score
  require(Seurat)
  de.genes <- FindMarkers(seurat_object, ident.1 = 1, ident.2 = 2, logfc.threshold = 1, min.pct = pct.1, recorrect_umi = FALSE)
  de.genes <- subset(de.genes, abs(de.genes$avg_log2FC) > 2)
  de.genes$p_val_adj[de.genes$p_val_adj == 0] <- 1e-20
  score <- de.genes$p_val_adj
  score <- -log10(score)
  score[score > 20] <- 20
  score <- sum(score)
  return(score)
}
