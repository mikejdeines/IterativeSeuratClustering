initial_clustering <- function(seurat_object, reduction.name = "pca", dims.use = 1:30){
  require(Seurat)
  require(scCustomize)
  seurat_object <- FindClusters(seurat_object, resolution = 0.01, algorithm = 4, method = "igraph", verbose = FALSE)
  objs <- SplitObject(seurat_object, "seurat_clusters")
  objs <- objs[order(names(objs))]
  samples <- lapply(objs, function(obj) leiden_clustering(obj, reduction.name = reduction.name, dims.use = dims.use))
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}

clustering_iteration <- function(seurat_object, min_score, cluster_size, pct.1, reduction.name = "pca", dims.use = 1:30){
  require(Seurat)
  require(scCustomize)
  Idents(seurat_object) <- seurat_object$leiden_clusters
  objs <- SplitObject(seurat_object, "leiden_clusters")
  objs <- objs[order(names(objs))]
  samples <- lapply(objs, function(obj) leiden_clustering(obj, score_limit = min_score, min_size = cluster_size, pct.1 = pct.1, reduction.name = reduction.name, dims.use = dims.use))
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}
